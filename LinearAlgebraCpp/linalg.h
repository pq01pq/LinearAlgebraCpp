#pragma once

#include "linalg_allocate.h"


#include <iostream>
#include <string>
#include <initializer_list>

namespace linalg {

	class Tensorr;
	class Matrixx;
	class Roww;
	class Vectorr;

	// Base class of vectors
	class Tensorr {
	public:
		virtual ~Tensorr() = default;
		virtual const size_t size() const = 0;

		virtual const std::string str() const = 0;
	protected:
		Tensorr() = default;
		static const double convertNegativeZero(const double value);
	};

	/*
	* When allocating, if sequence becomes larger than container size, it just pass by remaing values.
	* 
	* Index reference operators allows negative index as in python.
	* Range : -size ~ (size - 1)
	* 
	* '&' is a horizontal append operator and '|' is a vertical append operator.
	* Priority : & > | (follows default priority of 'and' and 'or')
	* 
	* Copy constructor, copy operator=, move constructor, move operator=, destructor : deprecated for Rule of Zero
	*/
	class Matrixx : public Tensorr, public Allocatablee {
		friend class Roww;
		friend class Vectorr;
	public:
		//Matrixx();
		Matrixx(const int height = 1, const int width = 1);
		Matrixx(const Matrixx& copyMatrix);
		//Matrixx(Matrixx&& moveMatrix) noexcept;
		explicit Matrixx(const Roww& copyRow);
		explicit Matrixx(const Vectorr& copyVector);
		virtual ~Matrixx() = default;
		void init(const int height = 1, const int width = 1); // throws std::length_error

		void reduce(); // == toEchelonForm + toReducedEchelonForm
		void toEchelonForm();
		void toReducedEchelonForm(); // throws std::logic_error

		bool isEchelonForm();

		Matrixx block(const size_t beginRow, const size_t beginCol,
			const size_t blockHeight, const size_t blockWidth) const; // throws std::out_of_range
		Matrixx inverse(); // throws std::logic_error, get inverse matrix of square matrix
		Matrixx transpose(const bool inplace = false);

		static Matrixx identity(const int length); // throws std::length_error, create elementary matrix(or unit matrix)
		//static Matrixx zero(const int height, const int width); // throws std::length_error, create zero matrix
		
		// Traditional array index reference method (only positive index)
		const Roww& operator[](const size_t row) const; // throws std::out_of_range
		Roww& operator[](const size_t row); // throws std::out_of_range

		// Modified index reference method (positive and negative index)
		const Roww& operator()(const int row) const; // throws std::out_of_range
		Roww& operator()(const int row); // throws std::out_of_range
		const double& operator()(const int row, const int col) const; // throws std::out_of_range
		double& operator()(const int row, const int col); // throws std::out_of_range
		// (recommended) operator() can catch both row and column index out of range exception.

		Allocatorr& operator<<(const double value); // 1st way to initialize entries
		Matrixx& operator=(const std::initializer_list<double> values); // 2nd way to initialize entries

		Matrixx operator+() const;
		Matrixx operator-() const;

		Matrixx& operator=(const Matrixx& rightMatrix);
		//Matrixx& operator=(Matrixx&& rightMatrix) noexcept;
		Matrixx& operator+=(const Matrixx& rightMatrix); // throws std::logic_error
		Matrixx& operator-=(const Matrixx& rightMatrix); // throws std::logic_error
		Matrixx& operator*=(const double multiplier);
		Matrixx& operator*=(const Matrixx& rightMatrix); // throws std::logic_error
		Matrixx& operator/=(const double divisor); // throws std::logic_error : divide by zero (for convenience)

		// Horizontal append operation
		Matrixx& operator&=(const Matrixx& rightMatrix); // throws std::logic_error
		Matrixx& operator&=(const Vectorr& rightVector); // throws std::logic_error
		// Vertical append operation
		Matrixx& operator|=(const Matrixx& lowerMatrix); // throws std::logic_error
		Matrixx& operator|=(const Roww& lowerRow); // throws std::logic_error

		// Vector equation operation
		friend Vectorr operator*(const Matrixx& leftMatrix, const Vectorr& rightVector); // throws std::logic_error

		friend bool operator==(const Matrixx& leftMatrix, const Matrixx& rightMatrix);

		virtual const size_t size() const;

		const size_t height() const;
		const size_t width() const;

		Roww getRow(const int row) const; // throws std::out_of_range
		Vectorr getColumn(const int col) const; // throws std::out_of_range

		virtual const std::string str() const override;

	protected:
		virtual void allocate(const size_t sequence, const double value) override;
	private:
		class Impl;
		std::unique_ptr<Impl> mMatrix;

		/*size_t mHeight, mWidth;
		std::vector<Roww> mRows;

		struct Pivot {
			size_t row, col;
			double entry;
		};*/

		/*
		* "컴퓨터 프로그램은 보통 한 열에서 가장 절댓값이 큰 성분을 추축으로 선정한다."
		* - [David C. Lay et al] Linear Algebra and Its Applications (선형대수학) 1.2 -
		*/
		//const Pivot findPivot(const size_t beginRow, const size_t beginCol) const; // Find largest absolute value of entries
		//const void replaceRowsUnder(const Pivot pivot); // Row replacing operation in forward phase
		//const Pivot getPivot(const size_t row) const; // Get existing pivot from row in echelon form matrix
		//const void replaceRowsOver(const Pivot pivot); // Row replacing operation in backward phase

		friend void swap(Matrixx& leftMatrix, Matrixx& rightMatrix) noexcept;
	};

	Matrixx operator+(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Matrixx operator-(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Matrixx operator*(const double multiplier, const Matrixx& rightMatrix);
	Matrixx operator*(const Matrixx& leftMatrix, const double multiplier);
	Matrixx operator*(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Matrixx operator/(const Matrixx& leftMatrix, const double divisor);

	// Horizontal append operation
	Matrixx operator&(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Matrixx operator&(const Matrixx& leftMatrix, const Vectorr& rightVector);
	Matrixx operator&(const Vectorr& leftVector, const Matrixx& rightMatrix);
	Matrixx operator&(const Vectorr& leftVector, const Vectorr& rightVector);
	// Vertical append operation
	Matrixx operator|(const Matrixx& upperMatrix, const Matrixx& lowerMatrix);
	Matrixx operator|(const Matrixx& upperMatrix, const Roww& lowerRow);
	Matrixx operator|(const Roww& upperRow, const Matrixx& lowerMatrix);
	Matrixx operator|(const Roww& upperRow, const Roww& lowerRow);

	bool operator!=(const Matrixx& leftMatrix, const Matrixx& rightMatrix);

	std::ostream& operator<<(std::ostream& outputStream, const Matrixx& outputMatrix);

	class Roww : public Tensorr, public Allocatablee {
		friend class Matrixx;
	public:
		//Roww();
		explicit Roww(const int size = 1);
		Roww(const Roww& copyRow);
		//Roww(Roww&& moveRow) noexcept;
		virtual ~Roww() = default;
		void init(const int size = 1);

		// Traditional array index reference method (only positive index)
		const double& operator[](const size_t col) const; // throws std::out_of_range
		double& operator[](const size_t col); // throws std::out_of_range

		// Modified index reference method (positive and negative index)
		const double& operator()(const int col) const; // throws std::out_of_range
		double& operator()(const int col); // throws std::out_of_range

		Allocatorr& operator<<(const double value); // 1st way to initialize entries
		Roww& operator=(const std::initializer_list<double> values); // 2nd way to initialize entries

		Roww operator+() const;
		Roww operator-() const;

		Roww& operator=(const Roww& rightRow);
		//Roww& operator=(Roww&& rightRow) noexcept;
		Roww& operator+=(const Roww& rightRow); // throws std::logic_error
		Roww& operator-=(const Roww& rightRow); // throws std::logic_error
		Roww& operator*=(const double multiplier);
		Roww& operator/=(const double divisor); // throws std::logic_error : divide by zero (for convenience)

		// Horizontal append operation
		Roww& operator&=(const Roww& rightRow);

		friend bool operator==(const Roww& leftRow, const Roww& rightRow);

		virtual const size_t size() const;
		const size_t width() const;

		virtual const std::string str() const override;
	protected:
		virtual void allocate(const size_t sequence, const double value) override;
	private:
		class Impl;
		std::unique_ptr<Impl> mRow;

		/*std::vector<double> mEntries;*/

		friend void swap(Roww& leftRow, Roww& rightRow) noexcept;
	};

	Roww operator+(const Roww& leftRow, const Roww& rightRow);
	Roww operator-(const Roww& leftRow, const Roww& rightRow);
	Roww operator*(const double multiplier, const Roww& rightRow);
	Roww operator*(const Roww& leftRow, const double multiplier);
	Roww operator/(const Roww& leftRow, const double divisor);

	// Horizontal append operation
	Roww operator&(const Roww& leftRow, const Roww& rightRow);

	bool operator!=(const Roww& leftRow, const Roww& rightRow);

	std::ostream& operator<<(std::ostream& outputStream, const Roww& outputRow);

	class Vectorr : public Tensorr, public Allocatablee {
		friend class Matrixx;
	public:
		//Vectorr();
		explicit Vectorr(const int size = 1);
		Vectorr(const Vectorr& copyVector);
		//Vectorr(Vectorr&& moveVector) noexcept;
		virtual ~Vectorr() = default;
		void init(const int size = 1);

		// Traditional array index reference method (only positive index)
		const double& operator[](const size_t row) const; // throws std::out_of_range
		double& operator[](const size_t row); // throws std::out_of_range

		// Modified index reference method (positive and negative index)
		const double& operator()(const int row) const; // throws std::out_of_range
		double& operator()(const int row); // throws std::out_of_range

		Allocatorr& operator<<(const double value); // 1st way to initialize entries
		Vectorr& operator=(const std::initializer_list<double> values); // 2nd way to initialize entries

		Vectorr operator+() const;
		Vectorr operator-() const;

		Vectorr& operator=(const Vectorr& rightVector);
		//Vectorr& operator=(Vectorr&& rightVector) noexcept;
		Vectorr& operator+=(const Vectorr& rightVector); // throws std::logic_error
		Vectorr& operator-=(const Vectorr& rightVector); // throws std::logic_error
		Vectorr& operator*=(const double multiplier);
		Vectorr& operator/=(const double divisor); // throws std::logic_error : divide by zero (for convenience)

		// Vertical append operation
		Vectorr& operator|=(const Vectorr& lowerVector);

		// Vector equation operation
		friend Vectorr operator*(const Matrixx& leftMatrix, const Vectorr& rightVector); // throws std::logic_error

		friend bool operator==(const Vectorr& leftVector, const Vectorr& rightVector);

		virtual const size_t size() const;
		const size_t height() const;

		virtual const std::string str() const override;
	protected:
		virtual void allocate(const size_t sequence, const double value) override;
	private:
		class Impl;
		std::unique_ptr<Impl> mVector;

		/*std::vector<double> mEntries;*/
		
		friend void swap(Vectorr& leftVector, Vectorr& rightVector) noexcept;
	};

	Vectorr operator+(const Vectorr& leftVector, const Vectorr& rightVector);
	Vectorr operator-(const Vectorr& leftVector, const Vectorr& rightVector);
	Vectorr operator*(const double multiplier, const Vectorr& rightVector);
	Vectorr operator*(const Vectorr& leftVector, const double multiplier);
	Vectorr operator/(const Vectorr& leftVector, const double divisor);

	// Vertical append operation
	Vectorr operator|(const Vectorr& upperVector, const Vectorr& lowerVector);

	bool operator!=(const Vectorr& leftVector, const Vectorr& rightVector);

	std::ostream& operator<<(std::ostream& outputStream, const Vectorr& outputVector);
}

#include "linalg_impl.h"