#ifndef LINSPACE_HPP
#define LINSPACE_HPP

#include <vector>

namespace HydrogenOrbitals
{
	/**
	* @brief linspace (vector of evenly spaced points).
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

		/// @brief Delta (step) value - the interval between consecutive points.
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
		* @param size  - the number of points between the start and end values.
		*
		* @exception std::invalid_argument - raised when size < 2.
		*/
		linspace(T start, T end, std::size_t size)
			: std::vector<T>()
			, delta(size > 1 ? (end - start) / (size - 1) : throw new std::invalid_argument("linspace size < 2."))
		{
			for (size_t i = 0; i < size; ++i)
			{
				this->push_back(start + i * delta);
			}
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

#endif /* !LINSPACE_HPP */