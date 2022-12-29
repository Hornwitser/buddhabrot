#include <algorithm>
#include <array>
#include <complex>
#include <iostream>
#include <random>
#include <string>
#include <string_view>
#include <vector>

#include <cstdlib>
#include <cstdint>

#include "lodepng.h"

#include "parse_args.hpp"
#include "srgb.hpp"


using std::int32_t;
using std::int64_t;

using PixelSize = int32_t;

struct Arguments {
    PixelSize width = 800;
    PixelSize height = 800;
    int64_t samples = 1'000'000;
    int64_t max_iterations = 1'000;

    std::string output_path = "render.png";
};

const std::array<OptionDescription<Arguments>, 6> option_descriptions = {
    "w", "width",          "width of the output image", &Arguments::width,
    "h", "height",         "height of the output image", &Arguments::height,
    "s", "samples",        "number of samples to compute", &Arguments::samples,
    "i", "max-iterations", "maximum number of iterations befor discarding path", &Arguments::max_iterations,
    "o", "output-path",    "path to output rendered PNG image", &Arguments::output_path,
    "?", "help",           "show this help", std::nullopt,
};

void plot_path(const Arguments& args, std::vector<int>& histogram, const std::vector<std::complex<float>>& path)
{
    for (const auto& point : path)
    {
        // XXX bad algorithm
        PixelSize x = (point.real() + 2.f) / 4.f * args.width;
        if (x < 0 || x >= args.width)
            continue;
        PixelSize y = (point.imag() + 2.f) / 4.f * args.height;
        if (y < 0 || y >= args.height)
            continue;
        histogram[x + y * args.width] += 1;
    }
}

int main(int argc, char* argv[])
{
    Arguments args;
    try
    {
        parse_args(argc, argv, args, option_descriptions);
    }
    catch (ParsingFailed)
    {
        return EXIT_FAILURE;
    }

    std::cout << "width: " << args.width << "\n";
    std::cout << "height: " << args.height << "\n";
    std::cout << "samples: " << args.samples << "\n";
    std::cout << "max_iterations: " << args.max_iterations << "\n";
    std::cout << "output_path: " << args.output_path << "\n";

    std::ranlux48_base engine;
    std::uniform_real_distribution<float> dist(-2, 2);
    std::vector<int> histogram(args.width * args.height);
    std::vector<std::complex<float>> path;
    for (int64_t i = 0; i < args.samples; i++)
    {
        path.clear();
        const std::complex<float> c(dist(engine), dist(engine));
        std::complex z = c;
        for (int64_t j = 0; j < args.max_iterations; j++)
        {
            path.push_back(z);
            z = z * z + c;
            if (std::norm(z) > 4)
            {
                plot_path(args, histogram, path);
                break;
            }
        }
    }


    int64_t total = std::reduce(histogram.begin(), histogram.end());
    std::cout << "total samples: " << total << std::endl;
    int64_t max = *std::max_element(histogram.begin(), histogram.end());
    std::cout << "lagest bucket size: " << max << std::endl;

    std::vector<uint8_t> image(args.width * args.height * 4);
    for (uint64_t i = 0; i < histogram.size(); i++)
    {
        uint8_t value = srgb_encoding_gamma((float)histogram[i] / max) * 255;
        // Write RGBA
        image[i * 4] = value;
        image[i * 4 + 1] = value;
        image[i * 4 + 2] = value;
        image[i * 4 + 3] = 255;
    }

    std::cout << "writing " << args.output_path << std::endl;
    unsigned error = lodepng::encode(args.output_path, image, args.width, args.height);
    if (error)
    {
        std::cerr << "encode error " << error << ": " << lodepng_error_text(error) << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
