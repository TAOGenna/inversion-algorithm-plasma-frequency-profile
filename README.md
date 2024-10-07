<div align="center">
  <img src="images/logo.jpg" alt="Plasma Frequency Profile" width="10%" />
  <h3>Inversion algorithm to retrieve the plasma frequency profile</h3>

  A quasi-parabolic approach for inverting ionograms
 
</div>
<p align="center">
  <img src="images/correct1.png" alt="Image 1" width="45%"/>
  <img src="images/correct2.png" alt="Image 2" width="45%"/>
</p>

## Overview
The inversion of ionograms into plasma frequency profiles is a longstanding problem, with several established approaches. However, most of these are implemented in closed software, often written in outdated programming languages. This project presents a Python implementation of an inversion algorithm based on quasi-parabolic layers. By alternating between these layers, we can construct a smooth plasma frequency profile while minimizing the least squares error between the generated and the original ionogram.

Refer to the accompanying PDF documentation for detailed information on the mathematics and algorithms used.

## How to use
This program expects `.SAO` files, which are commonly used for sharing remote sensing data. These files should be placed in the `/sao_files/` directory. To run the program, simply execute `python3 main.py` in your terminal. As the algorithm solves for the E and F layers, you will see a progress bar. The `avoid_date_list.txt` file lists dates with poor-quality ionograms that should be avoided; otherwise, the program will stop. The final image comparing the results should resemble the examples shown above.
<p align="center">
  <img src="images/progress_bar.png" alt="Image 3" width="95%"/>
</p>

## References
*The first paper has a lot of errors, but the underlying idea is correct. 
- L. Niu, L. Wen, C. Zhou, and M. Deng, ”A profile inversion method for vertical ionograms,” AIP Advances, vol. 14, no.6, p. 065034, Jun. 2024. doi: `10.1063/5.0208687`.
- J.E. Titheridge, ”A new method for the analysis of ionospheric h’(ƒ) records,” Journal of Atmospheric and Terrestrial Physics, vol. 21, no. 1, pp. 1-12, 1961. doi: `10.1016/0021-9169(61)90185-4`.
