This is the repository for barbi-map, a cumulative extinction map derived from the data provided by Barbillon et al. (2025, in prep.).

    ├── docs
    │   ├── build
    ├── README.md
    └── source
        ├── calc_ext.py
        ├── cumulative_ext_vaex.py
        └── show_somestars.py

There are three important files:
- `cumulative_ext_vaex.py` generates the extinction and reddening maps
- `calc_ext.py` contains the function to be called to recover the reddening and extinction from the map provided by `cumulative_ext_vaex.py`
- `show_somestars.py` performs comparison between **inferred** reddening in l22, bayestar, and this map

If you want to generate your own map, you will need `Yllari_selected_data.fits`. from the private folder data. If you want an out-of-the-box map, you can 
use `cumul_red_ag_barbillon.fits`.

There is also documentation compiled in `docs/build/html/`
To navigate it, simply download all the repo to your computer and open `docs/build/html/index.html` in any browser.