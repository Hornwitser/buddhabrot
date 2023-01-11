#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <iostream>
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
using std::uint32_t;
using std::int64_t;

using Clock = std::chrono::steady_clock;
using Microseconds = std::chrono::duration<double, std::micro>;
using Seconds = std::chrono::duration<double>;

struct Arguments {
    int64_t width = 800;
    int64_t height = 800;
    std::complex<float> centre = {0.f, 0.f};
    float zoom = 0.25f;
    int64_t max_iterations = 1'000;
    float samples_per_pixel = 1.f;
    std::optional<int64_t> samples = std::nullopt;
    std::optional<BoundingBox> sample_area = std::nullopt;
    int64_t mask_size = 1000;
    int32_t mask_edge_points = 4;
    std::optional<int64_t> mask_min_samples = std::nullopt;
    std::string mask_output_path = "";
    std::string point_density_output_path = "";

    bool escape_boundnary = false;
    std::optional<BoundingBox> output_area = std::nullopt;
    std::string output_path = "render.png";
};

const std::array<OptionDescription<Arguments>, 17> option_descriptions = {
    "w", "width",             "width of the output image", &Arguments::width,
    "h", "height",            "height of the output image", &Arguments::height,
    "c", "centre",            "coordinate of the centre of output image", &Arguments::centre,
    "z", "zoom",              "zoom level of output image", &Arguments::zoom,
    "S", "samples-per-pixel", "samples to compute normalized to output pixels", &Arguments::samples_per_pixel,
    "s", "samples",           "number of samples to compute (overrides -S)", &Arguments::samples,
    "a", "sample-area",       "area to make random samples in", &Arguments::sample_area,
    "i", "max-iterations",    "maximum number of iterations befor discarding path", &Arguments::max_iterations,
    "m", "mask-size",         "size of mandelbrot mask to use (0 for off)", &Arguments::mask_size,
    "e", "mask-edge-points",  "number of points to check along edges of the mask", &Arguments::mask_edge_points,
    "I", "mask-min-samples",  "minimum number of samples to check before bailing", &Arguments::mask_min_samples,
    "M", "mask-output-path",  "path to output mask for debugging", &Arguments::mask_output_path,
    "P", "point-output-path", "path to output point density map for debugging", &Arguments::point_density_output_path,
    "E", "escape",            "compute escape boundnary where magnitute only increases", &Arguments::escape_boundnary,
    "A", "output-area",       "area to show in output image (overrides -c, -z)", &Arguments::output_area,
    "o", "output-path",       "path to output rendered PNG image", &Arguments::output_path,
    "?", "help",              "show this help", std::nullopt,
};

struct Performance {
    int64_t mask_samples = 0;
    Clock::duration mask_time = Clock::duration(0);
    int64_t samples_input = 0;
    int64_t samples_mask = 0;
    int64_t samples_iterate = 0;
    int64_t samples_output = 0;
    int64_t points_input = 0;
    int64_t points_output = 0;
    Clock::duration buddhabrot_time = Clock::duration(0);
    Clock::duration compute_time = Clock::duration(0);
};

Mat<3, 3> area_to_image(const BoundingBox& area, int64_t width, int64_t height)
{
    return
        scale(width, -height) * translate(0.f, -1.f) *
        inverse(translate(area.min_x, area.min_y) * scale(area.max_x - area.min_x, area.max_y - area.min_y))
    ;
}

Mat<3, 3> image_to_area(int64_t width, int64_t height, const BoundingBox& area, bool centre = false)
{
    float shift = centre ? -0.5f : 0.f;
    return
        translate(area.min_x, area.min_y) * scale(area.max_x - area.min_x, area.max_y - area.min_y) *
        inverse(translate(shift, shift) * scale(width, -height) * translate(0.f, -1.f))
    ;
}

bool plot_path(
    const Arguments& args,
    const std::vector<std::complex<float>>& path,
    const Mat<3, 3>& transform,
    std::vector<uint32_t>& histogram,
    Performance& perf
) {
    int64_t points_output = 0;
    for (const auto& point : path)
    {
        ColVec<3> p = transform * ColVec<3>{point.real(), point.imag(), 1.f};
        if (p.x() < 0.f || args.width <= p.x() || p.y() < 0.f || args.height <= p.y())
            continue;
        histogram[(int64_t)(p.x()) + (int64_t)(p.y()) * args.width] += 1;
        points_output++;
    }

    perf.points_input += path.size();
    perf.points_output += points_output;
    return points_output;
}

void escape_boundnary(const Arguments& args, std::vector<uint32_t>& histogram, Performance& perf)
{
    Mat<3, 3> transform = image_to_area(args.width, args.height, *args.output_area, true);
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
    float d = std::sqrt(maximum_norm_distance(*args.output_area));
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

bool inside_mandelbrot(std::complex<float> c, int64_t max_iterations, float norm_limit)
{
    if (inside_cardioid(c) || inside_period2_bulb(c))
        return true;

    std::complex z = c;
    for (int64_t j = 0; j < max_iterations; j++)
    {
        z = z * z + c;
        if (std::norm(z) > norm_limit)
            return false;
    }

    return true;
}

void mandelbrot_mask(const Arguments& args, std::vector<bool>& mask, Performance& perf)
{
    Mat<3, 3> transform = image_to_area(args.mask_size, args.mask_size, *args.sample_area);
    //std::cout << "transform: " << transform << std::endl;
    float norm_limit = maximum_norm_distance(*args.output_area);

    static const std::array<ColVec<3>, 4> corner_offset = {
        0.f, 0.f, 0.f,
        1.f, 0.f, 0.f,
        1.f, 1.f, 0.f,
        0.f, 1.f, 0.f,
    };
    static const std::array<ColVec<3>, 4> segment_direction = {
        1.f, 0.f, 0.f,
        0.f, 1.f, 0.f,
        -1.f, 0.f, 0.f,
        0.f, -1.f, 0.f,
    };

    int64_t samples = 0;
    for (decltype(mask.size()) i = 0; i < mask.size(); i++)
    {
        ColVec<3> p = {(float)(i % args.mask_size), (float)(i / args.mask_size), 1.f};

        for (int segment = 0; segment < args.mask_edge_points; segment++)
            for (int corner = 0; corner < 4; corner++)
            {
                ColVec<3> p_edge = transform * (
                    p + corner_offset[corner] + (float)segment / args.mask_edge_points * segment_direction[corner]
                );
                std::complex c(p_edge.x(), p_edge.y());
                samples++;
                if (!inside_mandelbrot(c, args.max_iterations, norm_limit))
                {
                    mask[i] = false;
                    goto outside;
                }
            }

        mask[i] = true;
outside:;
    }

    perf.mask_samples += samples;
}

void buddhabrot(
    const Arguments& args,
    const std::vector<bool>& mask,
    std::vector<uint32_t>& histogram,
    std::optional<std::vector<uint32_t>>& point_density,
    Performance& perf
) {
    std::vector<std::complex<float>> path;

    Mat<3, 3> transform = area_to_image(*args.output_area, args.width, args.height);
    Mat<3, 3> mask_transform = image_to_area(args.mask_size, args.mask_size, *args.sample_area);
    //std::cout << "transform: " << transform << std::endl;

    // Using Roberts' sequence to select sample points.
    // See http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    const float g = 1.32471795724474602596f;
    const float a1 = 1.f / g;
    const float a2 = 1.f / (g * g);

    int64_t samples_per_mask_box = std::ceil((float)*args.samples / mask.size());
    std::cout << "samples per mask box: " << samples_per_mask_box << std::endl;
    float norm_limit = maximum_norm_distance(*args.output_area);
    for (decltype(mask.size()) i = 0; i < mask.size(); i++)
    {
        if (mask[i])
            continue;

        float x = i % args.mask_size;
        float y = i / args.mask_size;
        float xoff = 0.5f;
        float yoff = 0.5f;

        int64_t start_points = perf.points_output;
        bool has_points = false;
        for (int64_t j = 0; j < samples_per_mask_box; j++)
        {
            if (!has_points && *args.mask_min_samples && j > *args.mask_min_samples)
                break;

            perf.samples_mask++;
            ColVec<3> p = mask_transform * ColVec<3>{x + xoff, y + yoff, 1.f};
            xoff = std::fmod(xoff + a1, 1.f);
            yoff = std::fmod(yoff + a2, 1.f);
            const std::complex<float> c(p.x(), p.y());

            if (inside_cardioid(c) || inside_period2_bulb(c))
                continue;

            perf.samples_iterate++;
            path.clear();
            std::complex z = c;
            for (int64_t j = 0; j < args.max_iterations; j++)
            {
                path.push_back(z);
                z = z * z + c;
                if (std::norm(z) > norm_limit)
                {
                    if (plot_path(args, path, transform, histogram, perf))
                        has_points = true;
                    perf.samples_output++;
                    break;
                }
            }
        }
        if (point_density)
            (*point_density)[i] = perf.points_output - start_points;
    }

    perf.samples_input += mask.size() * samples_per_mask_box;
}

float area(const BoundingBox& box) {
    return (box.max_x - box.min_x) * (box.max_y - box.min_y);
}

template <typename T>
bool write_image(int64_t width, int64_t height, const std::vector<T>& histogram, const std::string& output_path)
{
    int64_t max = *std::max_element(histogram.begin(), histogram.end());
    std::cout << "lagest bucket size: " << max << std::endl;

    std::vector<uint8_t> image(width * height * 2);
    for (decltype(histogram.size()) i = 0; i < histogram.size(); i++)
    {
        // Note: This assumes passing 1.0f to srgb_encoding_gamma results in a value less than 1.0f being output.
        uint16_t value = srgb_encoding_gamma((float)histogram[i] / max) * 0x10000;
        image[i * 2] = value >> 8;
        image[i * 2 + 1] = value & 0xff;
    }

    std::cout << "writing " << output_path << std::endl;
    unsigned error = lodepng::encode(output_path, image, width, height, LCT_GREY, 16);
    if (error)
    {
        std::cerr << "encode error " << error << ": " << lodepng_error_text(error) << std::endl;
        return false;
    }

    return true;
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

    if (!args.output_area)
    {
        float half_width = 0.5f / args.zoom * (args.width > args.height ? (float)args.width / args.height : 1.f);
        float half_height = 0.5f / args.zoom * (args.width < args.height ? (float)args.height / args.width : 1.f);
        args.output_area = {
            args.centre.real() - half_width, args.centre.imag() - half_height,
            args.centre.real() + half_width, args.centre.imag() + half_height
        };
    }
    else
    {
        const BoundingBox& area = *args.output_area;
        args.zoom = 1.f / std::min(area.max_x - area.min_x, area.max_y - area.min_y);
        args.centre = {(area.min_x + area.max_x) / 2.f, (area.min_y + area.max_y) / 2.f};
    }

    if (!args.sample_area)
        args.sample_area = maximum_sample_area(args);

    if (!args.samples)
        args.samples = std::round(
            args.samples_per_pixel * args.height * args.width *
            area(*args.sample_area) / area(*args.output_area)
        );
    else
        args.samples_per_pixel =
            (float)*args.samples / args.height / args.width *
            area(*args.output_area) / area(*args.sample_area)
        ;

    if (!args.mask_min_samples)
        args.mask_min_samples = std::sqrt(*args.samples / std::max(args.mask_size * args.mask_size, (int64_t)1));

    std::cout << "width: " << args.width << "\n";
    std::cout << "height: " << args.height << "\n";
    std::cout << "centre: " << args.centre << "\n";
    std::cout << "zoom: " << args.zoom << "\n";
    std::cout << "samples_per_pixel: " << args.samples_per_pixel << "\n";
    std::cout << "samples: " << args.samples << "\n";
    std::cout << "sample_area: " << args.sample_area << "\n";
    std::cout << "mask_size: " << args.mask_size << "\n";
    std::cout << "mask_edge_points: " << args.mask_edge_points << "\n";
    std::cout << "mask_min_samples: " << args.mask_min_samples << "\n";
    std::cout << "max_iterations: " << args.max_iterations << "\n";
    std::cout << "escape_boundnary: " << args.escape_boundnary << "\n";
    std::cout << "output_area: " << args.output_area << "\n";
    std::cout << "output_path: " << args.output_path << "\n";

    std::vector<uint32_t> histogram(args.width * args.height);
    Performance perf;
    Clock::time_point start = Clock::now();
    if (args.escape_boundnary)
        escape_boundnary(args, histogram, perf);
    else
    {
        if (args.mask_size < 1)
            args.mask_size = 1;
        std::vector<bool> mask(args.mask_size * args.mask_size);
        if (args.mask_size > 1)
        {
            mandelbrot_mask(args, mask, perf);
            perf.mask_time = Clock::now() - start;
            std::cout << "mask_samples " << perf.mask_samples << "\n";
            std::cout << "mask_time: " << Seconds(perf.mask_time) << std::endl;
            std::cout << "mask megasamples/s: " << perf.mask_samples / Microseconds(perf.mask_time).count() << "\n";
            if (args.mask_output_path.size())
                write_image(args.mask_size, args.mask_size, mask, args.mask_output_path);
        }

        std::optional<std::vector<uint32_t>> point_density;
        if (args.point_density_output_path.size())
            point_density = std::vector<uint32_t>(args.mask_size * args.mask_size);

        Clock::time_point buddhabrot_start = Clock::now();
        buddhabrot(args, mask, histogram, point_density, perf);
        perf.buddhabrot_time = Clock::now() - buddhabrot_start;

        if (args.point_density_output_path.size())
            write_image(args.mask_size, args.mask_size, *point_density, args.point_density_output_path);
    }
    perf.compute_time = Clock::now() - start;

    std::cout << "samples_input: " << perf.samples_input << "\n";
    std::cout << "points_output: " << perf.points_output << "\n";
    std::cout << "samples_mask/input ratio: " << (float)perf.samples_mask / perf.samples_input << "\n";
    std::cout << "samples_iterate/mask ratio: " << (float)perf.samples_iterate / perf.samples_mask << "\n";
    std::cout << "samples_output/iterate ratio: " << (float)perf.samples_output / perf.samples_iterate << "\n";
    std::cout << "points_output/input ratio: " << (float)perf.points_output / perf.points_input << "\n";
    std::cout << "buddhabrot_time: " << Seconds(perf.buddhabrot_time) << "\n";
    std::cout << "input megasamples/s: " << perf.samples_input / Microseconds(perf.buddhabrot_time).count() << "\n";
    std::cout << "output megapoints/s: " << perf.points_output / Microseconds(perf.buddhabrot_time).count() << "\n";
    std::cout << "compute_time: " << Seconds(perf.compute_time) << "\n";
    std::cout << "buddhabrot/compute ratio: " << Seconds(perf.buddhabrot_time) / Seconds(perf.compute_time) << "\n";

    if (!write_image(args.width, args.height, histogram, args.output_path))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
