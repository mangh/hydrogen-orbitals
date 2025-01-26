#ifndef LINSPACE_HPP
#define LINSPACE_HPP

// #include "Numeral.hpp"
#include <vector>

namespace HydrogenOrbitals
{
    /**
    * @brief linspace vector
    * 
    */
    template<typename T>
    struct linspace : std::vector<T>
    {
        static_assert(std::is_arithmetic_v<T> == true);

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Properties
        ///

        // /// @brief Delta value.
        // const T delta;

        ///////////////////////////////////////////////////////////////////
        ///
        ///     Constructor(s)
        ///

        /**
        * @brief No default constructor!
        */
        linspace() = delete;

        // /**
        // * @brief linspace constructor.
        // * 
        // * @param start  - coordinate start value,
        // * @param end    - coordinate end value,
        // * @param count  - number of points (values) of the coordinate.
        // *
        // */
        // linspace(T start, T end, std::size_t count)
        //     : std::vector<T>(count)
        // {
        //     if (start >= end)
        //         throw new std::invalid_argument("linspace: start >= end.");

        //     if (count <= 0)
        //         throw new std::invalid_argument("linspace: count <= 0.");

        //     T delta = (end - start) / count;
        //     delta = (end + delta - start) / count;
            
        //     // for(size_t i = 0; i < this->size(); ++i) (*this)[i] = (start + i * delta);
        //     T x = start;
        //     for(auto &v : *this)
        //     {
        //         v = x;
        //         x += delta;
        //     }
        // }

        /**
        * @brief linspace constructor.
        * 
        * @param start  - coordinate start value,
        * @param end    - coordinate end value,
        * @param count  - number of points (values) of the coordinate.
        *
        */
        linspace(T start, T end, std::size_t count, bool end_included = true)
            : std::vector<T>()
        {
            if (start >= end)
                throw new std::invalid_argument("linspace: start >= end.");

            if (count <= 0)
                throw new std::invalid_argument("linspace: count <= 0.");

            T delta = (end - start) / count;
            T stop = end_included ? end + delta : end;
            while(start < stop)
            {
                this->push_back(start);
                start += delta;
            }
        }
    };
}

#endif /* !LINSPACE_HPP */