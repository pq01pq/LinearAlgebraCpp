#pragma once

#include "linalg.h"
#include "linalg_exception.h"

#include <vector>
#include <sstream>

namespace linalg {
	// Function implemented classes
	class Tensorr::Impl;
	class Matrixx::Impl;
	class Roww::Impl;
	class Vectorr::Impl;

	class Tensorr::Impl {
	public:
		virtual ~Impl() = default;

		virtual const size_t size() const = 0;

		virtual const std::string str() const = 0;
	protected:
		Impl() = default;

		static const double convertNegativeZero(const double value);
	};

	class Matrixx::Impl : public Tensorr::Impl {
		friend class Roww::Impl;
		friend class Vectorr::Impl;
	public:
		Impl(const int height = 1, const int width = 1);
		Impl(const Roww::Impl& copyRowImpl);
		Impl(const Vectorr::Impl& copyVectorImpl);
		virtual ~Impl() = default;
		void init(const int height = 1, const int width = 1); // throws std::length_error

		void reduce(); // == toEchelonForm + toReducedEchelonForm
		void toEchelonForm();
		void toReducedEchelonForm(); // throws std::logic_error

		bool isEchelonForm();

		Impl block(const size_t beginRow, const size_t beginCol,
			const size_t blockHeight, const size_t blockWidth) const; // throws std::out_of_range
		Impl inverse(); // throws std::logic_error, get inverse matrix of square matrix
		Impl transpose() const;

		static Impl identity(const int length); // throws std::length_error, create elementary matrix(or unit matrix)

		// Traditional array index reference method (only positive index)
		const Roww& operator[](const size_t row) const; // throws std::out_of_range
		Roww& operator[](const size_t row); // throws std::out_of_range

		// Modified index reference method (positive and negative index)
		const Roww& operator()(const int row) const; // throws std::out_of_range
		Roww& operator()(const int row); // throws std::out_of_range

		const double& operator()(const int row, const int col) const; // throws std::out_of_range
		double& operator()(const int row, const int col); // throws std::out_of_range
		// (recommended) operator() can catch both row and column index out of range exception.

		void allocate(const size_t sequence, const double value);
		void allocate(const std::initializer_list<double> values);

		Impl operator+() const;
		Impl operator-() const;

		Impl& operator+=(const Impl& rightMatrixImpl); // throws std::logic_error
		Impl& operator-=(const Impl& rightMatrixImpl); // throws std::logic_error
		Impl& operator*=(const double multiplier);
		Impl& operator*=(const Impl& rightMatrixImpl); // throws std::logic_error
		Impl& operator/=(const double divisor); // throws std::logic_error : divide by zero (for convenience)

		// Horizontal append operation
		Impl& operator&=(const Impl& rightMatrixImpl); // throws std::logic_error
		Impl& operator&=(const Vectorr::Impl& rightVectorImpl); // throws std::logic_error
		// Vertical append operation
		Impl& operator|=(const Impl& lowerMatrixImpl); // throws std::logic_error
		Impl& operator|=(const Roww::Impl& lowerRowImpl); // throws std::logic_error

		Impl operator+(const Impl& rightMatrixImpl) const;
		Impl operator-(const Impl& rightMatrixImpl) const;
		Impl operator*(const double multiplier) const;
		Impl operator*(const Impl& rightMatrixImpl) const;
		Impl operator/(const double divisor) const;

		// Vector equation operation
		Vectorr::Impl operator*(const Vectorr::Impl& rightVectorImpl) const; // throws std::logic_error

		// Horizontal append operation
		Impl operator&(const Impl& rightMatrixImpl) const;
		Impl operator&(const Vectorr::Impl& rightVectorImpl) const;
		
		// Vertical append operation
		Impl operator|(const Impl& lowerMatrixImpl) const;
		Impl operator|(const Roww::Impl& lowerRowImpl) const;
		
		bool operator==(const Impl& rightMatrixImpl) const;
		bool operator!=(const Impl& rightMatrixImpl) const;

		const Roww::Impl getRow(const int row) const;
		const Vectorr::Impl getColumn(const int col) const;

		const size_t height() const;
		const size_t width() const;
		virtual const size_t size() const override;

		virtual const std::string str() const override;
	private:
		struct Pivot {
			size_t row, col;
			double entry;
		};

		/*
		* "컴퓨터 프로그램은 보통 한 열에서 가장 절댓값이 큰 성분을 추축으로 선정한다."
		* - [David C. Lay et al] Linear Algebra and Its Applications (선형대수학) 1.2 -
		*/
		const Pivot findPivot(const size_t beginRow, const size_t beginCol) const; // Find largest absolute value of entries
		const void replaceRowsUnder(const Pivot pivot); // Row replacing operation in forward phase
		const Pivot getPivot(const size_t row) const; // Get existing pivot from row in echelon form matrix
		const void replaceRowsOver(const Pivot pivot); // Row replacing operation in backward phase

		//static Matrixx matrix(const Impl& matrixImpl);

		void swap(Impl& rightMatrixImpl) noexcept;

		size_t mHeight, mWidth;
		std::vector<Roww> mRows;
	};

	class Roww::Impl : public Tensorr::Impl {
		friend class Matrixx::Impl;
	public:
		Impl(const int size = 1);
		//Roww(Roww&& moveRow) noexcept;
		virtual ~Impl() = default;
		void init(const int size = 1);

		// Traditional array index reference
		const double& operator[](const size_t col) const; // throws std::out_of_range
		double& operator[](const size_t col); // throws std::out_of_range

		// Modified index reference method (positive and negative index)
		const double& operator()(const int col) const; // throws std::out_of_range
		double& operator()(const int col); // throws std::out_of_range

		void allocate(const size_t sequence, const double value);
		void allocate(const std::initializer_list<double> values);

		Impl operator+() const;
		Impl operator-() const;

		Impl& operator+=(const Impl& rightRow); // throws std::logic_error
		Impl& operator-=(const Impl& rightRow); // throws std::logic_error
		Impl& operator*=(const double multiplier);
		Impl& operator/=(const double divisor); // throws std::logic_error : divide by zero (for convenience)

		// Horizontal append operation
		Impl& operator&=(const Impl& rightRow);

		Impl operator+(const Impl& rightRowImpl) const;
		Impl operator-(const Impl& rightRowImpl) const;
		Impl operator*(const double multiplier) const;
		Impl operator/(const double divisor) const;

		// Horizontal append operation
		Impl operator&(const Impl& rightRowImpl) const;

		// Vertical append operation
		Matrixx::Impl operator|(const Matrixx::Impl& lowerMatrixImpl) const;
		Matrixx::Impl operator|(const Roww::Impl& lowerRowImpl) const;

		bool operator==(const Impl& rightRowImpl) const;
		bool operator!=(const Impl& rightRowImpl) const;

		const size_t width() const;
		virtual const size_t size() const override;

		virtual const std::string str() const override;
	private:
		void swap(Impl& rightRowImpl) noexcept;

		size_t mWidth;
		std::vector<double> mEntries;
	};

	

	

	class Vectorr::Impl : public Tensorr::Impl {
		friend class Matrixx::Impl;
	public:
		Impl(const int size = 1);
		//Vectorr(Vectorr&& moveVector) noexcept;
		virtual ~Impl() = default;
		void init(const int size = 1);

		// Traditional array index reference method (only positive index)
		const double& operator[](const size_t row) const; // throws std::out_of_range
		double& operator[](const size_t row); // throws std::out_of_range

		// Modified index reference method (positive and negative index)
		const double& operator()(const int row) const; // throws std::out_of_range
		double& operator()(const int row); // throws std::out_of_range

		void allocate(const size_t sequence, const double value);
		void allocate(const std::initializer_list<double> values);

		Impl operator+() const;
		Impl operator-() const;

		Impl& operator+=(const Impl& rightVectorImpl); // throws std::logic_error
		Impl& operator-=(const Impl& rightVectorImpl); // throws std::logic_error
		Impl& operator*=(const double multiplier);
		Impl& operator/=(const double divisor); // throws std::logic_error : divide by zero (for convenience)

		// Vertical append operation
		Impl& operator|=(const Impl& lowerVectorImpl);

		Impl operator+(const Impl& rightVectorImpl) const;
		Impl operator-(const Impl& rightVectorImpl) const;
		Impl operator*(const double multiplier) const;
		Impl operator/(const double divisor) const;

		// Horizontal append operation
		Matrixx::Impl operator&(const Matrixx::Impl& rightMatrixImpl) const;
		Matrixx::Impl operator&(const Impl& rightVectorImpl) const;

		// Vertical append operation
		Impl operator|(const Impl& lowerVectorImpl) const;

		bool operator==(const Impl& rightVectorImpl) const;
		bool operator!=(const Impl& rightVectorImpl) const;

		const size_t height() const;
		virtual const size_t size() const override;

		virtual const std::string str() const override;
	private:
		void swap(Impl& rightVectorImpl) noexcept;

		size_t mHeight;
		std::vector<double> mEntries;
	};
}