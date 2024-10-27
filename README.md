<div align="center">
  <img src="images/logo.jpg" alt="Plasma Frequency Profile" width="10%" />
  <h3>Inversion algorithm to retrieve the plasma frequency profile</h3>

  A quasi-parabolic approach for inverting ionograms
 
</div>

## Abstract
Inversion of ionograms from vertical measurements is significant for studying the ionospheric structure and ionospheric wave propagation, and it has attracted widespread attention. A model inversion algorithm based on the multivariate QP model (single quasi-parabolic ionospheric profile model) is proposed, which is different from the traditional QPS model (where the multiple quasi-parabolic segment ionospheric profile model uses quasi-parabolic or anti-parabolic models to represent E, valley, and F1 and F2 layers). In other words, the electron density height profile of each layer in the ionosphere is no longer described by a single QP model. Still, it is based on QP as the basic unit and characterized by a combination of multiple QP units, and the entire vertical ionospheric profile consists of a series of QP unit models. Moreover, in the case of the multivariate QP model, determining the parameters of each layer becomes more complex. Based on this model, the inversion of vertical ionograms was achieved. We basically reproduced the Niu et al. paper _A profile inversion method for vertical ionograms_, but offering a working-opensource code and improved upon the mathematics used since in the original paper the results seem a bit dubious. 

Refer to the accompanying <a href="Documentation.pdf">PDF documentation</a> for detailed information on the mathematics and algorithms used.
## Results
We used data provided by the Jicamarca Radio Observatory located in Lima, Peru. The data is a .SAO file thus providing all the information of the iongram such as the E and F layer critical frequency and the o-mode trace. We find that the program works for ionograms with E and F layer or daytime ionograms.

<p align="center">
  <img src="images/90.png" alt="Image 1" width="300px"/>
  <img src="images/92.png" alt="Image 2" width="300px"/>

  <img src="images/111.png" alt="Image 1" width="300px"/>
  <img src="images/130.png" alt="Image 2" width="300px"/>
</p>


## How to use
This program expects `.SAO` files, which are commonly used for sharing remote sensing data. These files should be placed in the `/sao_files/` directory. To run the program, simply execute `python3 main.py` in your terminal. As the algorithm solves for the E and F layers, you will see a progress bar. The `avoid_date_list.txt` file lists dates with poor-quality ionograms that should be avoided; otherwise, the program will stop. The final image comparing the results should resemble the examples shown above.
<p align="center">
  <img src="images/progress_bar.png" alt="Image 3" width="95%"/>
</p>

## References
- L. Niu, L. Wen, C. Zhou, and M. Deng, ”A profile inversion method for vertical ionograms,” AIP Advances, vol. 14, no.6, p. 065034, Jun. 2024. doi: `10.1063/5.0208687`.
- J.E. Titheridge, ”A new method for the analysis of ionospheric h’(ƒ) records,” Journal of Atmospheric and Terrestrial Physics, vol. 21, no. 1, pp. 1-12, 1961. doi: `10.1016/0021-9169(61)90185-4`.
