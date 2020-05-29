# ManifoldPlus: A Robust and Scalable Watertight Manifold Surface Generation Method for Triangle Soups
Advanced version of my previous Manifold algorithm from this [**repo**](https://github.com/hjwdzh/Manifold).

![Plane Fitting Results](https://github.com/hjwdzh/ManifoldPlus/raw/master/res/manifold-teaser.jpg)

### Dependencies
1. Eigen
2. LibIGL

### Installing prerequisites
```
git submodule update --init --recursive
```

### Quick examples
```
sh compile.sh
sh examples.sh
```

### Build
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
```

### Run
The input is a random triangle mesh in obj/off format. The output is a watertight manifold mesh in obj format.
```
./ManifoldPlus --input input.obj --output output.obj --depth 8
```
An example script is provided so that you can try several provided models. We convert inputs in data folder to outputs in results folder.

## Author
- [Jingwei Huang](mailto:jingweih@stanford.edu)

&copy; 2020 Jingwei Huang All Rights Reserved

**IMPORTANT**: If you use this code please cite the following (to provide) in any resulting publication:
```
@article{huang2020manifoldplus,
  title={ManifoldPlus: A Robust and Scalable Watertight Manifold Surface Generation Method for Triangle Soups},
  author={Huang, Jingwei and Zhou, Yichao and Guibas, Leonidas},
  journal={arXiv preprint arXiv:2005.11621},
  year={2020}
}
```
