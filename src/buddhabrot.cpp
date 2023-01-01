#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <iostream>
#include <random>
#include <string>
#include <string_view>
#include <vector>

#include <cstdlib>
#include <cstdint>

#include "lodepng.h"

#include "matrix.hpp"
#include "parse_args.hpp"
#include "srgb.hpp"


using std::int32_t;
using std::int64_t;

using PixelSize = int32_t;
using Clock = std::chrono::steady_clock;
using Microseconds = std::chrono::duration<double, std::micro>;
using Seconds = std::chrono::duration<double>;

struct Arguments {
    PixelSize width = 800;
    PixelSize height = 800;
    int64_t samples = 1'000'000;
    int64_t max_iterations = 1'000;
    std::optional<BoundingBox> sample_area = std::nullopt;

    bool escape_boundnary = false;
    BoundingBox output_area = { -2.f, -2.f, 2.f, 2.f };
    std::string output_path = "render.png";
};

const std::array<OptionDescription<Arguments>, 9> option_descriptions = {
    "w", "width",          "width of the output image", &Arguments::width,
    "h", "height",         "height of the output image", &Arguments::height,
    "s", "samples",        "number of samples to compute", &Arguments::samples,
    "a", "sample-area",    "area to make random samples in", &Arguments::sample_area,
    "i", "max-iterations", "maximum number of iterations befor discarding path", &Arguments::max_iterations,
    "E", "escape",         "compute escape boundnary where magnitute only increases", &Arguments::escape_boundnary,
    "A", "output-area",    "area to show in output image", &Arguments::output_area,
    "o", "output-path",    "path to output rendered PNG image", &Arguments::output_path,
    "?", "help",           "show this help", std::nullopt,
};

struct Performance {
    int64_t samples_input = 0;
    int64_t samples_iterate = 0;
    int64_t samples_output = 0;
    int64_t points_input = 0;
    int64_t points_output = 0;
    Clock::duration compute_time = Clock::duration(0);
};

void plot_path(
    const Arguments& args,
    const std::vector<std::complex<float>>& path,
    const Mat<3, 3>& transform,
    std::vector<int>& histogram,
    Performance& perf
) {
    int64_t points_output = 0;
    for (const auto& point : path)
    {
        ColVec<3> p = transform * ColVec<3>{point.real(), point.imag(), 1.f};
        if (p.x() < 0.f || args.width <= p.x() || p.y() < 0.f || args.height <= p.y())
            continue;
        histogram[(int)(p.x()) + (int)(p.y()) * args.width] += 1;
        points_output++;
    }

    perf.points_input += path.size();
    perf.points_output += points_output;
}

void escape_boundnary(const Arguments& args, std::vector<int>& histogram, Performance& perf)
{
    const BoundingBox& area = args.output_area;
    Mat<3, 3> transform =
        translate(area.min_x, area.min_y) * scale(area.max_x - area.min_x, area.max_y - area.min_y) *
        inverse(translate(-0.5f, -0.5f) * scale(args.width, args.height))
    ;
    //std::cout << "transform: " << transform << std::endl;

    for (decltype(histogram.size()) i = 0; i < histogram.size(); i++)
    {
        ColVec<3> p = transform * ColVec<3>{(float)(i % args.width), (float)(i / args.width), 1};

        std::complex c(p.x(), p.y());
        std::complex z = c;
        float norm_c = std::norm(c);
        for (int64_t j = 0; j < args.max_iterations; j++)
        {
            z = z * z + c;
            float norm_z = std::norm(z);
            if (!std::isfinite(norm_z))
                break;
            if (norm_z < norm_c)
            {
                histogram[i] = args.max_iterations - j;
                break;
            }
        }

        if (!std::isfinite(std::norm(z)) || std::norm(z) >= norm_c)
            histogram[i] = 0;
    }

    perf.samples_input += histogram.size();
}

/**
Calculates an upper bound on the squared distance from the origin point that
still includes all paths going through the bounding box.
*/
float maximum_norm_distance(BoundingBox box)
{
    // The structures around -2 + 0i can still reach anywhere in the fractal,
    // but outside a distance of 2 all paths rapidly increase away from the
    // origin. Thus a simple upper bound can be calculated by taking the larger
    // of the squared distance to the furthest point in the bounding and the
    // squared distance of 2.

    float x = std::max(std::abs(box.min_x), std::abs(box.max_x));
    float y = std::max(std::abs(box.min_y), std::abs(box.max_y));
    return std::max(2.f * 2.f, x * x + y * y);
}

/**
Calculates an upper bound for the size of the sample area needed to cover all
all points in the output area.
*/
BoundingBox maximum_sample_area(const Arguments& args)
{
    float d = std::sqrt(maximum_norm_distance(args.output_area));
    return {-d, -d, d, d};
}

bool inside_cardioid(std::complex<float> c)
{
    float p = std::sqrt(std::pow(c.real() - 1.f / 4.f, 2.f) + c.imag() * c.imag());
    return c.real() <= p - 2.f * p * p + 1.f / 4.f;
}

bool inside_period2_bulb(std::complex<float> c)
{
    return std::pow(c.real() + 1.f, 2.f) + std::pow(c.imag(), 2) <= 1.f / 16.f;
}

void buddhabrot(const Arguments& args, std::vector<int>& histogram, Performance& perf)
{
    std::ranlux48_base engine;
    std::uniform_real_distribution<float> x_dist(args.sample_area->min_x, args.sample_area->max_x);
    std::uniform_real_distribution<float> y_dist(args.sample_area->min_y, args.sample_area->max_y);
    std::vector<std::complex<float>> path;

    const BoundingBox& area = args.output_area;
    Mat<3, 3> transform =
        scale(args.width, args.height) *
        inverse(translate(area.min_x, area.min_y) * scale(area.max_x - area.min_x, area.max_y - area.min_y))
    ;
    //std::cout << "transform: " << transform << std::endl;

    float norm_limit = maximum_norm_distance(args.output_area);
    for (int64_t i = 0; i < args.samples; i++)
    {
        path.clear();
        const std::complex<float> c(x_dist(engine), y_dist(engine));

        if (inside_cardioid(c) || inside_period2_bulb(c))
            continue;
        perf.samples_iterate++;

        std::complex z = c;
        for (int64_t j = 0; j < args.max_iterations; j++)
        {
            path.push_back(z);
            z = z * z + c;
            if (std::norm(z) > norm_limit)
            {
                plot_path(args, path, transform, histogram, perf);
                perf.samples_output++;
                break;
            }
        }
    }

    perf.samples_input += args.samples;
}

int main(int argc, char* argv[])
{
    std::cout.setf(std::ios::boolalpha);
    std::cerr.setf(std::ios::boolalpha);

    Arguments args;
    try
    {
        parse_args(argc, argv, args, option_descriptions);
    }
    catch (ParsingFailed)
    {
        return EXIT_FAILURE;
    }

    if (!args.sample_area)
        args.sample_area = maximum_sample_area(args);

    std::cout << "width: " << args.width << "\n";
    std::cout << "height: " << args.height << "\n";
    std::cout << "samples: " << args.samples << "\n";
    std::cout << "sample_area: " << args.sample_area << "\n";
    std::cout << "max_iterations: " << args.max_iterations << "\n";
    std::cout << "escape_boundnary: " << args.escape_boundnary << "\n";
    std::cout << "output_area: " << args.output_area << "\n";
    std::cout << "output_path: " << args.output_path << "\n";

    std::vector<int> histogram(args.width * args.height);
    Performance perf;
    Clock::time_point start = Clock::now();
    if (args.escape_boundnary)
        escape_boundnary(args, histogram, perf);
    else
        buddhabrot(args, histogram, perf);
    perf.compute_time = Clock::now() - start;

    std::cout << "samples_input: " << perf.samples_input << "\n";
    std::cout << "points_output: " << perf.points_output << "\n";
    std::cout << "samples_iterate/input ratio: " << (float)perf.samples_iterate / perf.samples_input << "\n";
    std::cout << "samples_output/iterate ratio: " << (float)perf.samples_output / perf.samples_iterate << "\n";
    std::cout << "points_output/input ratio: " << (float)perf.points_output / perf.points_input << "\n";
    std::cout << "compute_time: " << Seconds(perf.compute_time) << "\n";
    std::cout << "input megasamples/s: " << perf.samples_input / Microseconds(perf.compute_time).count() << "\n";
    std::cout << "output megapoints/s: " << perf.points_output / Microseconds(perf.compute_time).count() << "\n";

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
