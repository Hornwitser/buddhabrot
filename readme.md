# Buddhabrot renderer in C++

A basic command line based multithreaded Buddhabrot fractal renderer that was written in C++ in 2 weeks.
Its main features are being resonably fast, able to plot about 200 million points per second on a quad core Intel Haswell CPU, and correctly handle zooming into the set without loosing points that should have ended up in the zoomed area.
The rendered result is automatically scaled and output as a 16-bit PNG with sRGB encoding gamma applied.

## Compiling

Currently only a modern Linux like environment is supported, you will need GCC 12 or newer to compile.
To build run `make` to compile and produce the `bin/buddhabrot` executable.

## Example Usage

```
buddhabrot --centre -1.2 0.706 --zoom 20 --samples-per-pixel 200 --output-path render.png
```

Renders a 20x zoomed in view centreed on the real coordinate -1.2 and imaginary coordinate 0.706 at a sample density of 200 samples per pixel in the output image and writes the result to render.png (overwriting if it already exists).
A zoom factor of 1 means render an image that's 1 unit wide in the real and imaginary axes, while 10 would be an image 0.1 units wide.
If the aspect ratio of the output is not square then the whichever is smaller of width and height is used as the reference.

Once you get the output render you will probably want to edit it with an image editor and at the very least adjust the mapping of pixel levels to something that makes better use of the full range as the default scaling produced by this tool is poor.

There are many options, see `buddhabrot --help` for a list of all of them.

## Future work

While I'm not planning on doing any more work on this project for the forseable future these are shortcommings I would like to have improved, ideas optimisations that could be done, and fetures I think would be cool:

### Better work thread scheduling

Right now work is partitioned out to threads based on an estimate of the "bogus operations per second" (bops for short) that can be computed by the system and scaled such that a work unit completes every half second.
This is problematic due to the wildly differeng amount of actual work that a bogus unit of work takes to compute and leads to a factor about 100 between the computation time near the edge of the sampling area where every sample escapes in one iteration and no samples hit the output and the middle of the sampling area where samples takes thousands of iterations and hit the output.
This in turn leads to bad thread usage becase one thread could end up with a giant work unit that takes several times longer to complete than the remaining work unitse divided on the rest of the threads, and so a bunch of cores stay idle at the end of the render.
It also leads to bad progress reporting as the progress is reported based on bogus work done rather than actual work done.

I don't have any good ideas for how to make this better, but making it better would improve the performance and stability of progress reporting.

### Better algorithm for ruling out starting locations

To make sure all sampling locations that could reach the output image are tested the current algorithm for selecting starting points for a zoomed in view is to indiscriminately sample small boxes from (-2, -2) to (2, 2) on the complex plane and give up early if no samples inside the box hit the output before a limit is reached.
This works well but leads to wasted work sampling locations like the box at (1.9, 1.9) to (2, 2) that would never hit the output.
It shouldn't be too hard to make a good filter estimate using few computations to rule out locations that definitevly don't hit the output.
This would also help improve the thread scheduling.

### Optimise for SIMD

The current code does not optimise well to SIMD instructions with GCC 12 which leads to less than ideal performance on modern CPUs.
I think speedups by as much as 2x could be acheived by making multiple calculations happen in parallel via AVX2 instructions.

### Store input parameters in the PNG

It would be useful to store all the parameters that was used to create a render as a text metadata block in the resulting PNG image, makes it much easier to track how a random late night render was created.

### Support clustered renderes over the network

This is just too cool to not have.
