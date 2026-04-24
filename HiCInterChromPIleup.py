import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cooltools.lib.numutils
import cooler
import pandas as pd
import os
import time
import argparse


# --- Helper Functions ---
def sum_normalization(matrix):
    total = matrix.sum()
    return matrix / total if total != 0 else matrix

def matrix_of_percentiles(matrix):
    matrix_sorted = np.sort(matrix, axis=None)
    percentiles = []
    for split in np.array_split(matrix_sorted, 100):
        percentiles.append(np.max(split))
        
    perc_df = pd.DataFrame({'percentile': np.linspace(1, 101, 100), 'value': percentiles})
    def return_percentile(x, df=perc_df):
        for el in df.value:
            if el >= x:
                return np.round(df[df.value == el].percentile.values[0] - 1)
    vfunc = np.vectorize(return_percentile, otypes=[float])
    return vfunc(matrix)

# --- Core Class ---
class RablPileup:
    def __init__(self, coolers, dim=50, skip_top=0, skip_bottom=0):
        self.dim = dim
        self.skip_top = skip_top
        self.skip_bottom = skip_bottom
        self.raw_matrix = self._compute_matrix(coolers)
        self.is_diff = False

    def _compute_matrix(self, coolers):
        total_matrix = np.zeros((self.dim, self.dim))
        for clr in coolers:
            matrix_sum = np.zeros((self.dim, self.dim))
            count = 0
            names = clr.chromnames
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    m = cooltools.lib.numutils.zoom_array(
                        clr.matrix(balance=False).fetch(names[i], names[j]),
                        same_sum=False, order=1, final_shape=(self.dim, self.dim)
                    )
                    matrix_sum += m
                    count += 1
            total_matrix += sum_normalization(matrix_sum / count)
        return total_matrix

    def cent_tel_score(self):
        cent = self.raw_matrix[0:int(self.dim/5), 0:int(self.dim/5)].sum()
        tel = self.raw_matrix[int(self.dim/5):, int(self.dim/2):].sum()
        return (cent + tel) / self.raw_matrix.sum() if self.raw_matrix.sum() != 0 else 0.0

    def _get_plot_matrix(self, log2=False, percentiles=False):
        mat = self.raw_matrix.copy()
        if log2:
            mat = np.log2(mat)
        mat = mat[self.skip_top:self.dim-self.skip_bottom, self.skip_top:self.dim-self.skip_bottom]
        if percentiles:
            mat = matrix_of_percentiles(mat)
        return mat

    def plot(self, plot_name=None, save=False, cmap='coolwarm', title=True,
            header_dict={'title': 'res: 100_000, size: 50', 'fontsize': 13},
            dpi=600, vmin=None, vmax=None, log2=False, percentiles=False,
            folder_to_save="./"):
        
        if self.is_diff:
            cmap = 'viridis'
        
        out_matrix = self._get_plot_matrix(log2, percentiles)
        if vmax is None:
            vmax = np.percentile(out_matrix, 99.5)
            
        print(f'out matrix stats:\n max value: {out_matrix.max()}\n min value: {out_matrix.min()}\n 99.5% percentile: {np.percentile(out_matrix, 99.5)}')

        fig, ax = plt.subplots()
        if title:
            plt.title(header_dict['title'], fontsize=header_dict['fontsize'])
        ax.set_aspect('equal')
        im = ax.matshow(out_matrix, cmap=cmap, vmin=vmin, vmax=vmax)
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        
        ticks = np.arange(0, self.dim-self.skip_top-self.skip_bottom, step=5)
        ax.set_xticks(ticks, [str(x) for x in ticks])
        ax.set_yticks(ticks, [str(x) for x in ticks])
        
        plt.colorbar(im, cax=cax)
        print(f'plot for {plot_name} was built!')
        
        if save:
            os.makedirs(folder_to_save, exist_ok=True)
            plt.savefig(os.path.join(folder_to_save, f"{plot_name}.png"), dpi=dpi, bbox_inches='tight')
            print('file was successfully saved to the hard drive')
        plt.show()
        time.sleep(1)

    def __sub__(self, other):
        if not isinstance(other, RablPileup):
            raise TypeError("Subtraction is only supported between two RablPileup objects.")
        new = RablPileup.__new__(RablPileup)
        new.dim = self.dim
        new.skip_top = self.skip_top
        new.skip_bottom = self.skip_bottom
        new.raw_matrix = self.raw_matrix - other.raw_matrix
        new.is_diff = True
        return new

# --- CLI Interface ---
def main():
    parser = argparse.ArgumentParser(description="Interchromosomal Pileup Analysis (Rabl Structure Visualization)")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # CREATE subcommand
    p_create = subparsers.add_parser("create", help="Create pileup from a list of cooler files")
    p_create.add_argument("coolers", nargs="+", help="List of .cool files")
    p_create.add_argument("--plot_name", default="pileup", help="Output plot name")
    p_create.add_argument("--save", action="store_true", help="Save plot to disk")
    p_create.add_argument("--log2", action="store_true", help="Apply log2 transform")
    p_create.add_argument("--percentiles", action="store_true", help="Convert to percentile ranks")
    p_create.add_argument("--dim", type=int, default=53, help="Matrix dimension")
    p_create.add_argument("--skip_top", type=int, default=2, help="Skip top rows/cols")
    p_create.add_argument("--skip_bottom", type=int, default=1, help="Skip bottom rows/cols")
    p_create.add_argument("--folder_to_save", default="./", help="Directory to save plots")

    # DIFF subcommand
    p_diff = subparsers.add_parser("diff", help="Compute difference between two sets of coolers")
    p_diff.add_argument("--coolers_a", nargs="+", required=True, help="First list of .cool files")
    p_diff.add_argument("--coolers_b", nargs="+", required=True, help="Second list of .cool files")
    p_diff.add_argument("--plot_name", default="diff_pileup", help="Output plot name")
    p_diff.add_argument("--save", action="store_true", help="Save plot to disk")
    p_diff.add_argument("--dim", type=int, default=53, help="Matrix dimension")
    p_diff.add_argument("--skip_top", type=int, default=2, help="Skip top rows/cols")
    p_diff.add_argument("--skip_bottom", type=int, default=1, help="Skip bottom rows/cols")
    p_diff.add_argument("--folder_to_save", default="./", help="Directory to save plots")

    args = parser.parse_args()

    if args.command == "create":
        coolers = [cooler.Cooler(f) for f in args.coolers]
        pileup = RablPileup(coolers, dim=args.dim, skip_top=args.skip_top, skip_bottom=args.skip_bottom)
        pileup.plot(plot_name=args.plot_name, save=args.save, log2=args.log2, percentiles=args.percentiles,
                    folder_to_save=args.folder_to_save)
                    
    elif args.command == "diff":
        coolers_a = [cooler.Cooler(f) for f in args.coolers_a]
        coolers_b = [cooler.Cooler(f) for f in args.coolers_b]
        
        p_a = RablPileup(coolers_a, dim=args.dim, skip_top=args.skip_top, skip_bottom=args.skip_bottom)
        p_b = RablPileup(coolers_b, dim=args.dim, skip_top=args.skip_top, skip_bottom=args.skip_bottom)
        
        diff_pileup = p_a - p_b
        
        out_mat = diff_pileup._get_plot_matrix()
        vmax = np.percentile(out_mat, 99.5)
        
        diff_pileup.plot(plot_name=args.plot_name, save=args.save, vmin=-vmax, vmax=vmax,
                        folder_to_save=args.folder_to_save)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()