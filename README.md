<div align="center">

  <h3>Inverse Algorithm to get the plasma frequency profile</h3>

  A quasi-parabolic approach for inverting ionograms.
 
</div>

## Overview
The problem of inverting ionograms into plasma frequency profiles is a longstanding problem with some known approaches. Most of them are closed software written on some antique language. This project is a Python implementation of an inversion algorithm based on quasi-parabolic layers. By alternating between these we are able the construct a smooth curve for the plasma frequency profile while fitting with an L2 error the produced vs the original ionogram.

## How to use
There are two parts in the `.ipynb` file where we have to put our data. First, at the beginning of the file, you have to fed the program with your ionogram data, that is the frequencies in Hertz and the virtual heights in meters of your ionogram. Finally, you have to provide where is the critical frequency of the E layer.

## Some results
<div align="center">
![Description of the image](correct1.png)
</div>

<div align="center">
![Description of the image](correct2.png)
</div>
