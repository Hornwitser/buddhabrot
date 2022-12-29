#include <algorithm>
#include <array>
#include <charconv>
#include <complex>
#include <concepts>
#include <iostream>
#include <random>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <cstdlib>
#include <cstdint>

#include "lodepng.h"

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

struct OptionDescription {
    std::string_view short_name;
    std::string_view long_name;
    std::string_view help;
    std::variant<std::monostate, int32_t Arguments::*, int64_t Arguments::*, std::string Arguments::*> field;
};

const std::array<OptionDescription, 6> option_descriptions = {
    "w", "width",          "width of the output image", &Arguments::width,
    "h", "height",         "height of the output image", &Arguments::height,
    "s", "samples",        "number of samples to compute", &Arguments::samples,
    "i", "max-iterations", "maximum number of iterations befor discarding path", &Arguments::max_iterations,
    "o", "output-path",    "path to output rendered PNG image", &Arguments::output_path,
    "?", "help",           "show this help", std::monostate(),
};

struct ParsingFailed {};

void show_help(char* argv[])
{
    std::cerr << "Usage: " << argv[0] << " [options]\n";
    Arguments default_arguments;
    for (const auto& opt : option_descriptions)
        std::cerr << "  -" << opt.short_name << ", --" << opt.long_name << " " << opt.help << "\n";
}

const OptionDescription& find_short_opt(std::string_view option)
{
    for (auto& description : option_descriptions)
        if (option == description.short_name)
            return description;
    std::cerr << "Invalid short option " << option << std::endl;
    throw ParsingFailed{};
}

const OptionDescription& find_long_opt(std::string_view option)
{
    for (auto& description : option_descriptions)
        if (option == description.long_name)
            return description;
    std::cerr << "Invalid long option " << option << std::endl;
    throw ParsingFailed{};
}

void parse_value(std::monostate, Arguments&, std::string_view)
{ }

template <std::integral T>
void parse_value(T Arguments::* field, Arguments& args, std::string_view value)
{
    std::from_chars_result result = std::from_chars(value.begin(), value.end(), args.*field);
    if (result.ec == std::errc::result_out_of_range)
    {
        std::cerr << "Failed to parse int, input " << value << " out of range" << std::endl;
        throw ParsingFailed{};
    }
    else if (result.ptr != value.end())
    {
        std::cerr << "Failed to parse int, input " << value << " is not an integer" << std::endl;
        throw ParsingFailed{};
    }
}

void parse_value(std::string Arguments::* field, Arguments& args, std::string_view value)
{
    args.*field = value;
}


int parse_opts(int argc, char* argv[], Arguments& args, int argument_pos)
{
    std::string_view option_argument = argv[argument_pos++];
    if (option_argument.length() < 2)
    {
        std::cerr << "Invalid argument " << option_argument << std::endl;
        throw ParsingFailed{};
    }

    const OptionDescription* description;
    if (option_argument[1] == '-')
    {
        std::string_view option = option_argument.substr(2);
        description = &find_long_opt(option);
    }
    else
    {
        std::string_view option = option_argument.substr(1, 2);
        description = &find_short_opt(option);
    }

    if (description->long_name == "help")
    {
        show_help(argv);
        throw ParsingFailed{};
    }

    std::string_view value;
    if (!std::holds_alternative<std::monostate>(description->field))
    {
        if (argument_pos >= argc)
        {
            std::cerr << "Expected argument containing value for " << description->long_name << std::endl;
            throw ParsingFailed{};
        }
        else
            value = argv[argument_pos];
    }

    std::visit([&](auto& field) { parse_value(field, args, value); }, description->field);
    return 1;
}

void parse_args(int argc, char* argv[], Arguments& args)
{
    int argument_pos = 1;
    while (argument_pos < argc)
    {
        std::string_view current_arg = argv[argument_pos++];
        if (current_arg.length() < 1)
            continue; // Ignore empty arguments
        if (current_arg[0] == '-')
            argument_pos += parse_opts(argc, argv, args, argument_pos - 1);
        else
        {
            // Unsupported positional argument
            std::cerr << "Unknown positional argument " << current_arg << std::endl;
            throw ParsingFailed{};
        }
    }
}

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
        parse_args(argc, argv, args);
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
        uint8_t value = histogram[i] * 255 / max;
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
