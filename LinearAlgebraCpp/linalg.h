#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <initializer_list>
#include "linalg_exception.h"

namespace linalg {

	class Allocator;
	class Allocatable;
	class LinalgContainer;
	class Matrixx;
	class Roww;
	class Vectorr;

	/*
	* Allocator class is used for allocating values into linear algebra containers.
	* It allocates values when it encounters operator '<<' and ',' with increating sequence.
	*/
	class Allocator {
		friend class Allocatable;
	public:
		Allocator(Allocatable& target, const int sequence);
		~Allocator();
		Allocator& operator,(const double value);
		
	private:
		Allocatable& target;
		int sequence;
	};

	// Mix-in class of linear algebra containers(== Java interface)
	class Allocatable {
		friend class Allocator;
	protected:
		inline virtual void allocate(const int sequence, const double value) {}
	};

	// Base class of linear algebra containers
	class LinalgContainer : protected Allocatable {
	public:
		~LinalgContainer();
		const int size() const;

		inline virtual const std::string str() const { return "NaN"; }
	protected:
		LinalgContainer();
		LinalgContainer(int size);

		int mSize;

		const double convertNegativeZero(const double value) const;
	};

	/*
	* When allocating, if sequence becomes larger than container size, it just pass by remaing values.
	* 
	* Index reference operators allows negative index as in python.
	* Range : -size ~ (size - 1)
	* 
	* 
	*/
	class Matrixx : public LinalgContainer {
		friend class Roww;
		friend class Vectorr;
	public:
		Matrixx();
		Matrixx(const int height, const int width);
		Matrixx(const Matrixx& copyMatrix);
		explicit Matrixx(const Roww& copyRow);
		explicit Matrixx(const Vectorr& copyVector);
		~Matrixx();
		void init(const int height, const int width); // throws std::length_error

		void reduce(); // == toEchelonForm + toReducedEchelonForm
		void toEchelonForm();
		void toReducedEchelonForm(); // throws std::logic_error

		friend bool isEchelonForm(const Matrixx& matrix);

		Matrixx block(const int beginRow, const int beginCol,
			const int blockHeight, const int blockWidth) const; // throws std::out_of_range

		static Matrixx identity(const int length); // throws std::length_error, create elementary matrix(or unit matrix)
		static Matrixx zero(const int height, const int width); // throws std::length_error, create zero matrix
		Matrixx inverse(); // throws std::logic_error, get inverse matrix of square matrix
		
		Roww& operator[](const int row); // throws std::out_of_range
		const Roww& operator[](const int row) const; // throws std::out_of_range
		double& operator()(const int row, const int col) const; // throws std::out_of_range
		// (recommended) operator() can catch both row and column index out of range exception.

		Matrixx& operator=(std::initializer_list<double> values); // 1st way to initialize entries
		Allocator& operator<<(const double value); // 2nd way to initialize entries

		Matrixx operator+() const;
		Matrixx operator-() const;

		Matrixx& operator=(const Matrixx& rightMatrix);
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

		virtual const int getHeight() const;
		virtual const int getWidth() const;

		Roww getRow(const int row) const; // throws std::out_of_range
		Vectorr getColumn(const int col) const; // throws std::out_of_range

		virtual const std::string str() const override;

	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		int mHeight, mWidth;
		std::unique_ptr<Roww[]> mRows;

		struct Pivot {
			int row, col;
			double entry;
		};

		const Pivot findPivot(const int beginRow, const int beginCol) const; // Find largest absolute value of entries
		const void replaceRowsUnder(const Pivot pivot); // Row replacing operation in forward phase
		const Pivot getPivot(const int row) const; // Get existing pivot from row in echelon form matrix
		const void replaceRowsOver(const Pivot pivot); // Row replacing operation in backward phase

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

	class Roww : public LinalgContainer {
		friend class Matrixx;
		friend class Vectorr;
	public:
		Roww();
		explicit Roww(const int size);
		Roww(const Roww& copyRow);
		~Roww();
		void init(const int size);

		double& operator[](const int col); // throws std::out_of_range
		const double& operator[](const int col) const; // throws std::out_of_range

		Roww& operator=(std::initializer_list<double> values); // 1st way to initialize entries
		Allocator& operator<<(const double value); // 2nd way to initialize entries

		Roww operator+() const;
		Roww operator-() const;

		Roww& operator=(const Roww& rightRow);
		Roww& operator+=(const Roww& rightRow); // throws std::logic_error
		Roww& operator-=(const Roww& rightRow); // throws std::logic_error
		Roww& operator*=(const double multiplier);
		Roww& operator/=(const double divisor); // throws std::logic_error : divide by zero (for convenience)

		// Horizontal append operation
		Roww& operator&=(const Roww& rightRow);

		friend bool operator==(const Roww& leftRow, const Roww& rightRow);

		virtual const std::string str() const override;
	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		std::unique_ptr<double[]> mEntries;

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

	class Vectorr : public LinalgContainer {
		friend class Matrixx;
		friend class Roww;
	public:
		Vectorr();
		explicit Vectorr(const int size);
		Vectorr(const Vectorr& copyVector);
		~Vectorr();
		void init(const int size);

		double& operator[](const int row); // throws std::out_of_range
		const double& operator[](const int row) const; // throws std::out_of_range

		Vectorr& operator=(std::initializer_list<double> values); // 1st way to initialize entries
		Allocator& operator<<(const double value); // 2nd way to initialize entries

		Vectorr operator+() const;
		Vectorr operator-() const;

		Vectorr& operator=(const Vectorr& rightVector);
		Vectorr& operator+=(const Vectorr& rightVector); // throws std::logic_error
		Vectorr& operator-=(const Vectorr& rightVector); // throws std::logic_error
		Vectorr& operator*=(const double multiplier);
		Vectorr& operator/=(const double divisor); // throws std::logic_error : divide by zero (for convenience)

		// Vertical append operation
		Vectorr& operator|=(const Vectorr& lowerVector);

		// Vector equation operation
		friend Vectorr operator*(const Matrixx& leftMatrix, const Vectorr& rightVector); // throws std::logic_error

		friend bool operator==(const Vectorr& leftVector, const Vectorr& rightVector);

		virtual const std::string str() const override;
	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		std::unique_ptr<double[]> mEntries;
		
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
