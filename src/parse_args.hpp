#include <charconv>
#include <concepts>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>


struct ParsingFailed {};

template <typename T>
std::ostream& operator << (std::ostream& oss, std::optional<T> opt)
{
    if (opt)
        oss << *opt;
    else
        oss << "auto";
    return oss;
}

struct BoundingBox {
    float min_x, min_y;
    float max_x, max_y;
};

std::ostream& operator << (std::ostream& oss, BoundingBox box)
{
    return oss << box.min_x << " " << box.min_y << " " << box.max_x << " " << box.max_y;
}

template <typename T>
struct OptionDescription {
    std::string_view short_name;
    std::string_view long_name;
    std::string_view help;
    std::optional<std::variant<
        bool T::*,
        int32_t T::*,
        int64_t T::*,
        std::string T::*,
        BoundingBox T::*,
        std::optional<BoundingBox> T::*
    >> field;
};

template <typename T, typename Opts>
void show_help(char* argv[], const Opts& option_descriptions)
{
    std::cerr << "Usage: " << argv[0] << " [options]\n";
    T default_arguments;
    for (const auto& opt : option_descriptions)
    {
        std::cerr << "  -" << opt.short_name << ", --" << opt.long_name << " " << opt.help;
        if (opt.field)
            std::visit([&](auto field) {
                std::cerr << " (default " << default_arguments.*field << ")";
            }, *opt.field);
        std::cerr << "\n";
    }
}

template <typename T, typename Opts>
const OptionDescription<T>& find_short_opt(std::string_view option, const Opts& option_descriptions)
{
    for (auto& description : option_descriptions)
        if (option == description.short_name)
            return description;
    std::cerr << "Invalid short option " << option << std::endl;
    throw ParsingFailed{};
}

template <typename T, typename Opts>
const OptionDescription<T>& find_long_opt(std::string_view option, const Opts& option_descriptions)
{
    for (auto& description : option_descriptions)
        if (option == description.long_name)
            return description;
    std::cerr << "Invalid long option " << option << std::endl;
    throw ParsingFailed{};
}

template <typename T>
void parse_number(std::string_view value, T& field)
{
    std::from_chars_result result = std::from_chars(value.begin(), value.end(), field);
    if (result.ec == std::errc::result_out_of_range)
    {
        std::cerr << "Failed to parse number, input " << value << " out of range" << std::endl;
        throw ParsingFailed{};
    }
    if (result.ptr != value.end())
    {
        std::cerr << "Failed to parse number, input " << value << " is not valid" << std::endl;
        throw ParsingFailed{};
    }
}

template <typename T>
int parse_value(int, char* [], T&, int, bool* field)
{
    *field = true;
    return 0;
}

template <std::integral Int, typename T>
int parse_value(int argc, char* argv[], T&, int argp, Int* field)
{
    if (argp >= argc)
    {
        std::cerr << "Expected integer but got end of argument list" << std::endl;
        throw ParsingFailed{};
    }
    parse_number(argv[argp], *field);
    return 1;
}

template <typename T>
int parse_value(int argc, char* argv[], T&, int argp, std::string* field)
{
    if (argp >= argc)
    {
        std::cerr << "Expected string but got end of argument list" << std::endl;
        throw ParsingFailed{};
    }
    *field = argv[argp];
    return 1;
}

template <typename T>
int parse_value(int argc, char* argv[], T&, int argp, BoundingBox* field)
{
    if (argp + 4 > argc)
    {
        std::cerr << "Expected 4 floats but got end of argument list" << std::endl;
        throw ParsingFailed{};
    }
    parse_number(argv[argp], (*field).min_x);
    parse_number(argv[argp+1], (*field).min_y);
    parse_number(argv[argp+2], (*field).max_x);
    parse_number(argv[argp+3], (*field).max_y);
    return 4;
}

template <typename T, typename U>
int parse_value(int argc, char* argv[], T& args, int argp, std::optional<U>* field)
{
    if (argp < argc && std::string_view(argv[argp]) == "auto")
    {
        field->reset();
        return 1;
    }
    field->emplace();
    return parse_value(argc, argv, args, argp, &(field->value()));
}


template <typename T, typename Opts>
int parse_opts(int argc, char* argv[], T& args, int argp, const Opts& option_descriptions)
{
    std::string_view option_argument = argv[argp++];
    if (option_argument.length() < 2)
    {
        std::cerr << "Invalid argument " << option_argument << std::endl;
        throw ParsingFailed{};
    }

    const OptionDescription<T>* description;
    if (option_argument[1] == '-')
    {
        std::string_view option = option_argument.substr(2);
        description = &find_long_opt<T>(option, option_descriptions);
    }
    else
    {
        std::string_view option = option_argument.substr(1, 2);
        description = &find_short_opt<T>(option, option_descriptions);
    }

    if (description->long_name == "help")
    {
        show_help<T>(argv, option_descriptions);
        throw ParsingFailed{};
    }

    if (!description->field)
        return 0;

    return std::visit(
        [&](auto& field) { return parse_value(argc, argv, args, argp, &(args.*field)); }, *description->field
    );
}

template <typename T, typename Opts>
void parse_args(int argc, char* argv[], T& args, const Opts& option_descriptions)
{
    int argp = 1;
    while (argp < argc)
    {
        std::string_view current_arg = argv[argp++];
        if (current_arg.length() < 1)
            continue; // Ignore empty arguments
        if (current_arg[0] == '-')
            argp += parse_opts(argc, argv, args, argp - 1, option_descriptions);
        else
        {
            // Unsupported positional argument
            std::cerr << "Unknown positional argument " << current_arg << std::endl;
            throw ParsingFailed{};
        }
    }
}
