# occam1dtem
Toy Occam inversion code for 1D TEM (central loop) method in Matlab
## DATA
The example data are from the forward modelling of a synthetic model (examples/true.mod)
The forward modelling algorithm is an approximate method called ABFM (Adaptive Born Forward Mapping) from Christensen, 2002:

Christensen N. B. (2002). A generic 1-D imaging method for transient electromagnetic data. Geophysics, 67(2), 438–447. doi:10.1190/1.1468603
see also my toy occam code for magnetotellurics:

[https://github.com/dong-hao/occam1dmt]

## USAGE
See example/testbench.m for a simple demonstration on how to load the data and call the inversion code. 

## something like a disclaimer

This was one of many toy codes I fiddled with when I was a student - I hope this could be useful to our students nowadays in the EM community. 
Those who want to try this script are free to use it on academic/educational cases. But of course, I cannot guarantee the script to be working properly and calculating correctly (although I wish so). Have you any questions or suggestions, please feel free to contact me (but don't you expect that I will reply quickly!).  

## HOW TO GET IT
```
git clone https://github.com/dong-hao/occam1dtem/ your_local_folder
```

## UNITS
The internal scale here is log10(Simens/m) - instead of linear scale for both conductivity and apparent conductivity. The layer depth is in (linear) metres, while the dB/dt unit is in A/m^2. 

## ERRORS    
Currently the internal error here is standard deviation.

## HOW TO GET UPDATED
```
cd to_you_local_folder
git pull 
```

## Contact

DONG Hao –  donghao@cugb.edu.cn

China University of Geosciences, Beijing 

Distributed under the GPL v3 license. See ``LICENSE`` for more information.

[https://github.com/dong-hao/occam1dtem]

## Contributing

Those who are willing to contribute are welcomed to try - but I probably won't have the time to review the commits frequently (not that I would expect there will be any). 

1. Fork it (<https://github.com/dong-hao/occam1dtem/fork>)
2. Create your feature branch (`git checkout -b feature/somename`)
3. Commit your changes (`git commit -am 'Add some features'`)
4. Push to the branch (`git push origin feature/somename`)
5. Create a new Pull Request - lather, rinse, repeat 
