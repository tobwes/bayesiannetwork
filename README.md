# bayesiannetwork

This is a fork of gneuralnetwork 0.8 with the intention to conduct a feasibility study of bayesian networks.

I have come quite far, at least it is running now and you can see some results in the RStatistics folder (basically, there is one plot displaying the results for the curve_fitting_quadratic_function test).

I have also tried the test "curve_fitting_square_root" and it runs well, too.


## How to compile

There is a new dependency, namely lapacke. To install this in debian, run

`apt-get install liblapacke liblapacke-dev`

To compile, run

```
git clone https://github.com/tobwes/bayesiannetwork.git
autoreconf --install
./configure
make
```

Then, to run some tests, e.g.

`src/gneural_network tests/curve_fitting_quadratic_function.input`

If you have installed R, the RStatistics folder has a saved workspace with some scripts to display the results, e.g. stats_curve_fitting_quadratic()

Please be informed that this code is, due to my poor knowlege of C, far from perfect. There are, for example, most definitely some memory leaks in the code.

## The future

Depending on the interest of the developers of gneural network in this project, I am willing to port it back into the gneural network code. If they are not interested, I will probably still patch the code and remove unnecessary dependencies, but most likely development will fade out slowly...