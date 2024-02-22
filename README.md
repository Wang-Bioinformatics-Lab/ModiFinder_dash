# ModiFinder dash

This repository provides codes for the web version of:
``` ModiFinder: Tandem Mass Spectral Alignment Enables Structural Modification Site Localization ```

_Mohammad Reza Zare Shahneh, Michael Strobel, Giovanni Andrea Vitale, Christian Geibel, Vanessa V Phelan, Daniel Petras, Allegra T Aron, Yasin El Abiead, Neha Garg, Mingxun Wang_

## Running ModiFinder Locally

1. Make sure the checkout this repository and checkout the submodule with the following commands

    ```
    git submodule init
    git submodule update
    ```
1. Add SIRIUS data (SIRIUS folder in tht main directory) - optional

1. run the following
    ```
    make server-compose-interactive
    ```

You should be able to access the server now on http://localhost:5999


