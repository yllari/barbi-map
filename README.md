This is the repository for barbi-map, a cumulative extinction map derived from the data provided by Barbillon et al. (2025, in prep.).

There are three important files:

    ```cumulative_ext_vaex.py``` generates the extinction and reddening maps
    ```calc_ext.py``` contains the function to be called to recover the reddening and extinction from the map provided by ```cumulative_ext_vaex.py```
    ```show_somestars.py``` performs comparison between **inferred** reddening in l22, bayestar, and this map

There is also documentation compiled in docs/build/
To navigate it, simply open docs/build/index.html in any browser.