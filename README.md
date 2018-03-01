# taylor-law-mortality

`taylor-law-mortality` provides the R source code used in the paper: 

_Cohen, Bohk-Ewald, Rau (2018). Gompertz, Makeham, and Siler models explain Taylor’s law in human mortality data. Demographic Research 38(29): 773-842_

## How to run

The code needs to be executed in four steps, as defined by `step-*.R` files in the root directory.

### Prerequisites

- The `.R` files in `./dumps/` contain empty lists due to copyrights. You will need to put in mortality data yourself. See comments in [`step-1-mortality-parameter-estimation.R`](./step-1-mortality-parameter-estimation.R).

- Make sure you have configured the correct paths (in each of the four `step-*.R` files).

### Execution

Run the `step-*.R` scripts in the root directory in the prescribed order.

## How to cite

If you use this code for academic research, please cite the above paper.

## How to contribute

Please note that this source code is an academic project. We welcome any issues and pull requests.

## License

The source code of `taylor-law-mortality` is published under the [GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.en.html). 
