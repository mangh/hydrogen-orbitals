#ifndef MANGH_LINSPACE_HPP
#define MANGH_LINSPACE_HPP

#include <vector>

namespace mangh::math
{
	/**
	* @brief linspace (a vector of evenly spaced numbers of type T).
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

		/// @brief Delta (step) value - the interval between consecutive numbers.
		const T delta;

		///////////////////////////////////////////////////////////////////
		///
		///     Constructor(s)
		///

		/// @brief No default constructor!
		linspace() = delete;

		/**
		* @brief linspace constructor.
		*
		* @param start - the start value,
		* @param end   - the end value,
		* @param size  - the number of points between the `start` and `end` values.
		*
		* @retval 
		*
		* @li - for `size > 1`: vector of `size` elements evenly spaced between the `start` and `end` values,
		*
		* @li - for `size = 1`: one-element vector (with the only element set to `start`),
		*
		* @li - for `size < 1`: empty vector.
		*
		* 
		*/
		linspace(T start, T end, std::size_t size)
			: std::vector<T>()
			, delta(size > 1 ? (end - start) / (size - 1) : T{})
		{
			for (size_t i = 0; i < size; ++i) this->push_back(start + i * delta);
		}

		/**
		* @brief linspace copy constructor.
		*
		* @param other - linspace object to copy from.
		*
		*/
		linspace(const linspace& other)
			: std::vector<T>(other)
			, delta(other.delta)
		{
		}
	};
}

#endif /* !MANGH_LINSPACE_HPP */