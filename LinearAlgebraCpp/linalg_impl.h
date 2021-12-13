#pragma once

#include "linalg.h"
#include "linalg_exception.h"

#include <vector>
#include <sstream>

namespace linalg {
	// Function implemented classes
	class Matrixx::Impl;
	class Roww::Impl;
	class Vectorr::Impl;

	class Matrixx::Impl {
	public:
		Impl(const int height = 1, const int width = 1);
		Impl(const Roww& copyRow);
		Impl(const Vectorr& copyVector);
		~Impl() = default;
		void init(const int height = 1, const int width = 1); // throws std::length_error

		void reduce(); // == toEchelonForm + toReducedEchelonForm
		void toEchelonForm();
		void toReducedEchelonForm(); // throws std::logic_error

		bool isEchelonForm();

		Matrixx block(const size_t beginRow, const size_t beginCol,
			const size_t blockHeight, const size_t blockWidth) const; // throws std::out_of_range
		Matrixx inverse(); // throws std::logic_error, get inverse matrix of square matrix
		Matrixx transpose() const;

		static Matrixx identity(const int length); // throws std::length_error, create elementary matrix(or unit matrix)

		void allocate(const size_t sequence, const double value);
		void allocate(const std::initializer_list<double> values);

		Matrixx negative() const;
		void innerAdd(const Matrixx& rightMatrix);
		void innerSub(const Matrixx& rightMatrix);
		void innerMul(const double multiplier);
		void innerMatMul(const Matrixx& rightMatrix);
		void innerDiv(const double divisor);
		void innerHorizontalAppend(const Matrixx& rightMatrix);
		void innerHorizontalAppend(const Vectorr& rightVector);
		void innerVerticalAppend(const Matrixx& lowerMatrix);
		void innerVerticalAppend(const Roww& lowerRow);

		Vectorr vectorEquation(const Vectorr& rightVector);

		bool equals(const Matrixx& rightMatrix) const;

		const Roww& row(const size_t row) const;
		Roww& row(const size_t row);
		const Roww& row(const int row) const;
		Roww& row(const int row);
		const double& get(const int row, const int col) const;
		double& get(const int row, const int col);

		const Roww getRow(const int row) const;
		const Vectorr getColumn(const int col) const;

		const size_t height() const;
		const size_t width() const;

		const std::string str() const;
	private:
		size_t mHeight, mWidth;
		std::vector<Roww> mRows;

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

		static Matrixx matrix(const Impl& matrixImpl);

		void swap(Impl& rightMatrixImpl) noexcept;
	};

	class Roww::Impl {
	public:
		Impl(const int size = 1);
		//Roww(Roww&& moveRow) noexcept;
		~Impl() = default;
		void init(const int size = 1);

		void allocate(const size_t sequence, const double value);
		void allocate(const std::initializer_list<double> values);

		Roww negative() const;
		void innerAdd(const Roww& rightRow);
		void innerSub(const Roww& rightRow);
		void innerMul(const double multiplier);
		void innerDiv(const double divisor);
		void innerHorizontalAppend(const Roww& rightRow);

		bool equals(const Roww& rightRow) const;

		const double& get(const size_t col) const;
		double& get(const size_t col);
		const double& get(const int col) const;
		double& get(const int col);

		const size_t width() const;

		const std::string str() const;
	private:
		size_t mWidth;
		std::vector<double> mEntries;

		void swap(Impl& rightRowImpl) noexcept;
	};

	class Vectorr::Impl {
	public:
		Impl(const int size = 1);
		//Vectorr(Vectorr&& moveVector) noexcept;
		~Impl() = default;
		void init(const int size = 1);

		void allocate(const size_t sequence, const double value);
		void allocate(const std::initializer_list<double> values);

		Vectorr negative() const;
		void innerAdd(const Vectorr& rightVector);
		void innerSub(const Vectorr& rightVector);
		void innerMul(const double multiplier);
		void innerDiv(const double divisor);
		void innerVerticalAppend(const Vectorr& lowerVector);

		bool equals(const Vectorr& rightVector) const;

		const double& get(const size_t row) const;
		double& get(const size_t row);
		const double& get(const int row) const;
		double& get(const int row);

		const size_t height() const;

		const std::string str() const;
	private:
		size_t mHeight;
		std::vector<double> mEntries;

		void swap(Impl& rightVectorImpl) noexcept;
	};
}