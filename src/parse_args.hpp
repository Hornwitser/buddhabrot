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
struct OptionDescription {
    std::string_view short_name;
    std::string_view long_name;
    std::string_view help;
    std::optional<std::variant<int32_t T::*, int64_t T::*, std::string T::*>> field;
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

template <std::integral Int, typename T>
void parse_value(Int T::* field, T& args, std::string_view value)
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

template <typename T>
void parse_value(std::string T::* field, T& args, std::string_view value)
{
    args.*field = value;
}


template <typename T, typename Opts>
int parse_opts(int argc, char* argv[], T& args, int argument_pos, const Opts& option_descriptions)
{
    std::string_view option_argument = argv[argument_pos++];
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

    if (argument_pos >= argc)
    {
        std::cerr << "Expected argument containing value for " << description->long_name << std::endl;
        throw ParsingFailed{};
    }

    std::visit([&](auto& field) { parse_value(field, args, argv[argument_pos]); }, *description->field);
    return 1;

}

template <typename T, typename Opts>
void parse_args(int argc, char* argv[], T& args, const Opts& option_descriptions)
{
    int argument_pos = 1;
    while (argument_pos < argc)
    {
        std::string_view current_arg = argv[argument_pos++];
        if (current_arg.length() < 1)
            continue; // Ignore empty arguments
        if (current_arg[0] == '-')
            argument_pos += parse_opts(argc, argv, args, argument_pos - 1, option_descriptions);
        else
        {
            // Unsupported positional argument
            std::cerr << "Unknown positional argument " << current_arg << std::endl;
            throw ParsingFailed{};
        }
    }
}
