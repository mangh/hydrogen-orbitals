#ifndef MANGH_NPY_HPP
#define MANGH_NPY_HPP

#include <complex>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "cppfeatures.hpp"

// The .npy files are written as little-endian ('<' in descr), straight from memory:
// the platform must be little-endian (checked at compile time).
#if CPP_VERSION >= 2020
#  include <bit>
   static_assert(std::endian::native == std::endian::little, "The .npy writer assumes a little-endian platform.");
#elif defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)    // GCC, Clang
   static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "The .npy writer assumes a little-endian platform.");
#elif defined(_WIN32)                                                 // MSVC: Windows is little-endian on all its architectures
#else
#  error "Cannot determine the byte order: the .npy writer assumes a little-endian platform."
#endif

namespace mangh { namespace npy
{

    /// @brief NumPy type descriptor ("descr") for the element type T.
    template<typename T>
    constexpr const char* descr()
    {
        if constexpr (std::is_same_v<T, double>)                    return "<f8";
        else if constexpr (std::is_same_v<T, std::complex<double>>) return "<c16";
        else if constexpr (std::is_same_v<T, std::int32_t>)         return "<i4";
        else static_assert(!sizeof(T), "Unsupported .npy element type.");
    }

    /**
     * @brief Writes the .npy header (format version 1.0) for an array of the element type T and the given shape.
     *
     * @remarks The data written afterwards must be in C (row-major) order.
     * See https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html
     */
    template<typename T>
    void write_header(std::ostream& os, std::initializer_list<std::size_t> shape)
    {
        std::string dims;
        for (std::size_t d : shape) dims += std::to_string(d) + ", ";
        // Python tuple syntax: "(n,)" for 1-D, "(a, b, c)" otherwise
        if (shape.size() == 1) dims.pop_back();
        else if (shape.size() > 1) dims.resize(dims.size() - 2);

        std::string header = std::string("{'descr': '") + descr<T>() + "', 'fortran_order': False, 'shape': (" + dims + "), }";

        // magic (6) + version (2) + header length (2) + header, padded with spaces and ended by '\n' to a multiple of 64 bytes
        constexpr std::size_t preamble = 10;
        header.append(63 - (preamble + header.size()) % 64, ' ');
        header += '\n';

        const auto len = static_cast<std::uint16_t>(header.size());
        os.write("\x93NUMPY\x01\x00", 8);
        os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        os.write(header.data(), static_cast<std::streamsize>(header.size()));
    }

    /**
     * @brief Writes a contiguous array (in C order) as a .npy file.
     *
     * @param file  - output file path,
     * @param data  - pointer to the first element,
     * @param shape - array shape (the product of the dimensions is the number of elements).
     *
     * @exception std::runtime_error - when the file cannot be opened or written.
     */
    template<typename T>
    void save(const std::filesystem::path& file, const T* data, std::initializer_list<std::size_t> shape)
    {
        std::size_t count = 1;
        for (std::size_t d : shape) count *= d;

        std::ofstream os(file, std::ios::binary);
        if (!os) throw std::runtime_error("Cannot open " + file.string() + " for writing");
        write_header<T>(os, shape);
        os.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(count * sizeof(T)));
        os.close();
        if (!os) throw std::runtime_error("Failed to write " + file.string());
    }
}}

#endif /* !MANGH_NPY_HPP */
