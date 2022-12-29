#include <cmath>

/**
Encodes a normalized linear light intensity in the range 0-1 into sRGB intensity
see https://en.wikipedia.org/wiki/SRGB#From_CIE_XYZ_to_sRGB
*/
inline float srgb_encoding_gamma(float linear_input)
{
    if (linear_input <= 0.0031308f)
    {
        return 12.92f * linear_input;
    }
    return 1.055f * std::pow(linear_input, 1.f / 2.4f) - 0.055f;
}
