# HiCInterChromPIleup

A Python tool for computing and visualizing interchromosomal Hi-C contact pileups, optimized for Rabl chromosome architecture analysis.

## Repository Files
| File | Description |
|------|-------------|
| `HiCInterChromPIleup.py` | Core script containing the `RablPileup` class and a CLI for generating pileup & difference plots. |
| `usage_example.ipynb` | Jupyter notebook with interactive examples, custom parameter tuning, and visualization workflows. |

## Installation & Dependencies
```bash
git clone https://github.com/tsoiegor/HiCInterChromPIleup.git
cd HiCInterChromPIleup
pip install numpy matplotlib cooltools cooler pandas
```

## CLI Usage
The script exposes two commands: `create` (single condition) and `diff` (difference between two conditions).

### Generate a Pileup Plot
```bash
python HiCInterChromPIleup.py create sample1.mcool::/resolutions/100000 sample2.mcool::/resolutions/100000 \
    --plot_name rabl_pileup \
    --save \
    --dim 53 \
    --skip_top 2 \
    --skip_bottom 1
```

### Generate a Difference Plot
```bash
python HiCInterChromPIleup.py diff \
    --coolers_a aux1.mcool::/resolutions/100000 aux2.mcool::/resolutions/100000 \
    --coolers_b ctrl1.mcool::/resolutions/100000 ctrl2.mcool::/resolutions/100000 \
    --plot_name aux_vs_ctrl_diff \
    --save \
    --dim 53 \
    --skip_top 2 \
    --skip_bottom 1
```

## Notes
- For `.mcool` files, always append `::/resolutions/<bin_size>` (e.g., `::/resolutions/100000`).
- All plotting options (`--log2`, `--percentiles`, `--cmap`, `--dpi`, etc.) are available. Run `python HiCInterChromPIleup.py create --help` or `diff --help` for the full argument list.
- See `usage_example.ipynb` for programmatic/interactive usage of the `RablPileup` class.
