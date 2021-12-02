#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "linalg_exception.h"

namespace linalg {

	class Allocator;
	class Allocatable;
	class Matrixx;
	class Roww;
	class Vectorr;

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

	class Allocatable {
		friend class Allocator;
	protected:
		virtual void allocate(const int sequence, const double value) {};

		const double convertNegativeZero(const double value) const;
	};

	class Matrixx : private Allocatable {
		friend class Roww;
		friend class Vectorr;
	public:
		Matrixx(const int height, const int width);
		Matrixx(const Matrixx& copyMatrix);
		explicit Matrixx(const Roww& copyRow);
		explicit Matrixx(const Vectorr& copyVector);
		~Matrixx();

		void reduce();
		void toEchelonForm();
		void toReducedEchelonForm();

		friend bool isEchelonForm(const Matrixx& matrix);

		Matrixx block(const int beginRow, const int beginCol, const int blockHeight, const int blockWidth) const;

		static Matrixx identity(const int length);
		static Matrixx zero(const int height, const int width);
		Matrixx inverse();

		Roww& operator[](const int row);
		const Roww& operator[](const int row) const;
		double& operator()(const int row, const int col) const;

		Allocator& operator<<(const double value);

		Matrixx operator+() const;
		Matrixx operator-() const;

		Matrixx& operator=(const Matrixx& rightMatrix);
		Matrixx& operator+=(const Matrixx& rightMatrix);
		Matrixx& operator-=(const Matrixx& rightMatrix);
		Matrixx& operator*=(const double multiplier);
		Matrixx& operator*=(const Matrixx& rightMatrix);

		Matrixx& operator&=(const Matrixx& rightMatrix);
		Matrixx& operator&=(const Vectorr& rightVector);
		Matrixx& operator|=(const Matrixx& lowerMatrix);
		Matrixx& operator|=(const Roww& lowerRow);

		friend Vectorr operator*(const Matrixx& leftMatrix, const Vectorr& rightVector);

		const int getHeight() const;
		const int getWidth() const;

		Roww getRow(const int row) const;
		Vectorr getColumn(const int col) const;

		const std::string str() const;

	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		Roww* rows;
		int height, width;

		struct Pivot {
			int row, col;
			double entry;
		};

		const Pivot findPivot(const int beginRow, const int beginCol) const;
		const void replaceRowsUnder(const Pivot pivot);
		const Pivot getPivot(const int row) const;
		const void replaceRowsOver(const Pivot pivot);

		friend void swap(Matrixx& leftMatrix, Matrixx& rightMatrix) noexcept;
	};

	Matrixx operator+(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Matrixx operator-(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Matrixx operator*(const double multiplier, const Matrixx& rightMatrix);
	Matrixx operator*(const Matrixx& leftMatrix, const double multiplier);
	Matrixx operator*(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Vectorr operator*(const Matrixx& leftMatrix, const Vectorr& rightVector);

	Matrixx operator&(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Matrixx operator&(const Matrixx& leftMatrix, const Vectorr& rightVector);
	Matrixx operator&(const Vectorr& leftVector, const Matrixx& rightMatrix);
	Matrixx operator&(const Vectorr& leftVector, const Vectorr& rightVector);

	Matrixx operator|(const Matrixx& upperMatrix, const Matrixx& lowerMatrix);
	Matrixx operator|(const Matrixx& upperMatrix, const Roww& lowerRow);
	Matrixx operator|(const Roww& upperRow, const Matrixx& lowerMatrix);
	Matrixx operator|(const Roww& upperRow, const Roww& lowerRow);

	bool operator==(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	bool operator!=(const Matrixx& leftMatrix, const Matrixx& rightMatrix);

	std::ostream& operator<<(std::ostream& outputStream, const Matrixx& outputMatrix);

	class Roww : private Allocatable {
		friend class Matrixx;
		friend class Vectorr;
	public:
		Roww() = default;
		explicit Roww(const int width);
		Roww(const Roww& copyRow);
		~Roww();

		double& operator[](const int col);
		const double& operator[](const int col) const;

		Allocator& operator<<(const double value);

		Roww operator+() const;
		Roww operator-() const;

		Roww& operator=(const Roww& rightRow);
		Roww& operator+=(const Roww& rightRow);
		Roww& operator-=(const Roww& rightRow);
		Roww& operator*=(const double multiplier);

		Roww& operator&=(const Roww& rightRow);

		const int getWidth() const;

		const std::string str() const;
	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		double* entries;
		int width;
		
		void init(const int width);

		friend void swap(Roww& leftRow, Roww& rightRow) noexcept;
	};

	Roww operator+(const Roww& leftRow, const Roww& rightRow);
	Roww operator-(const Roww& leftRow, const Roww& rightRow);
	Roww operator*(const double multiplier, const Roww& rightRow);
	Roww operator*(const Roww& leftRow, const double multiplier);

	Roww operator&(const Roww& leftRow, const Roww& rightRow);

	bool operator==(const Roww& leftRow, const Roww& rightRow);
	bool operator!=(const Roww& leftRow, const Roww& rightRow);

	std::ostream& operator<<(std::ostream& outputStream, const Roww& outputRow);

	class Vectorr : private Allocatable {
		friend class Matrixx;
		friend class Roww;
	public:
		Vectorr() = default;
		explicit Vectorr(const int height);
		Vectorr(const Vectorr& copyVector);
		~Vectorr();

		double& operator[](const int row);
		const double& operator[](const int row) const;

		Allocator& operator<<(const double value);

		Vectorr operator+() const;
		Vectorr operator-() const;

		Vectorr& operator=(const Vectorr& rightVector);
		Vectorr& operator+=(const Vectorr& rightVector);
		Vectorr& operator-=(const Vectorr& rightVector);
		Vectorr& operator*=(const double multiplier);

		Vectorr& operator|=(const Vectorr& lowerVector);

		friend Vectorr operator*(const Matrixx& leftMatrix, const Vectorr& rightVector);

		const int getHeight() const;

		const std::string str() const;
	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		double* entries;
		int height;

		void init(const int height);

		friend void swap(Vectorr& leftVector, Vectorr& rightVector) noexcept;
	};

	Vectorr operator+(const Vectorr& leftVector, const Vectorr& rightVector);
	Vectorr operator-(const Vectorr& leftVector, const Vectorr& rightVector);
	Vectorr operator*(const double multiplier, const Vectorr& rightVector);
	Vectorr operator*(const Vectorr& leftVector, const double multiplier);
	Vectorr operator*(const Matrixx& leftMatrix, const Vectorr& rightVector);

	Vectorr operator|(const Vectorr& upperVector, const Vectorr& lowerVector);

	bool operator==(const Vectorr& leftVector, const Vectorr& rightVector);
	bool operator!=(const Vectorr& leftVector, const Vectorr& rightVector);

	std::ostream& operator<<(std::ostream& outputStream, const Vectorr& outputVector);
}
