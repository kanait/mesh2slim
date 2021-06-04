# mesh2slim

mesh2slim is a research code for hierarchical approximation of SLIM surfaces from polygonal meshes, which includes an implementation of the following paper:

Takashi Kanai, Yutaka Ohtake, Kiwamu Kase: “Hierarchical Error-Driven Approximation of Implicit Surfaces from Polygonal Meshes”, Proc. 4th Eurographics/ACM SIGGRAPH Symposium on Geometry Processing, pp.21-30, 2006.

This software was originally developed in 2005-2006 and was renovated in 2021 so as to build successfully on CMake environment.

## Getting Started

This software is a command-line based application. First you execute "git clone" with with "--recursive" option.

```
git clone https://github.com/kanait/mesh2slim.git --recursive
```

You then execute in mesh2slim directory:

```
% cd mesh2slim
% mkdir build
% cd build
% cmake ..
% make
```
Then an executable "mesh2slim" is created if compilation completed successfully.

## How to Use

"mesh2slim" is used to convert from .obj (Wavefront OBJ) file to .slim2t (SLIM) file.
```
% mesh2slim in.obj out.slim2t
```

To display SLIM (.slim2t) file, use [slimviewGPU](https://github.com/kanait/slimviewGPU).

## Prerequisites

The following libraries are required for successfully compiling this software.

### [vecmath-cpp](https://github.com/yuki12/vecmath-cpp)
### [Eigen](https://gitlab.com/libeigen/eigen)

## Authors

* **[Takashi Kanai](https://graphics.c.u-tokyo.ac.jp/hp/en/)** - The University of Tokyo

## License

This software is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
