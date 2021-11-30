#pragma once
#include <iostream>
#include <string>

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

		Matrixx block(const int beginRow, const int beginCol, const int blockHeight, const int blockWidth) const;

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
			double value;
		};

		friend void swap(Matrixx& leftMatrix, Matrixx& rightMatrix) noexcept;
	};

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

		const int getWidth() const;

		const std::string str() const;
	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		double* entries;
		int width;
		
		void init(int width);

		friend void swap(Roww& leftRow, Roww& rightRow) noexcept;
	};

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

	Matrixx identityMatrix(const int length);
	Matrixx zeroMatrix(const int height, const int width);

	Matrixx operator+(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Matrixx operator-(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Matrixx operator*(const double multiplier, const Matrixx& rightMatrix);
	Matrixx operator*(const Matrixx& leftMatrix, const double multiplier);
	Matrixx operator*(const Matrixx& leftMatrix, const Matrixx& rightMatrix);

	Matrixx operator&(const Matrixx& leftMatrix, const Matrixx& rightMatrix);
	Matrixx operator&(const Matrixx& leftMatrix, const Vectorr& rightVector);
	Matrixx operator&(const Vectorr& leftVector, const Matrixx& rightMatrix);
	Matrixx operator&(const Vectorr& leftVector, const Vectorr& rightVector);

	Matrixx operator|(const Matrixx& upperMatrix, const Matrixx& lowerMatrix);
	Matrixx operator|(const Matrixx& upperMatrix, const Roww& lowerRow);
	Matrixx operator|(const Roww& upperRow, const Matrixx& lowerMatrix);
	Matrixx operator|(const Roww& upperRow, const Roww& lowerRow);

	Roww operator+(const Roww& leftRow, const Roww& rightRow);
	Roww operator-(const Roww& leftRow, const Roww& rightRow);
	Roww operator*(const double multiplier, const Roww& rightRow);
	Roww operator*(const Roww& leftRow, const double multiplier);

	Vectorr operator+(const Vectorr& leftVector, const Vectorr& rightVector);
	Vectorr operator-(const Vectorr& leftVector, const Vectorr& rightVector);
	Vectorr operator*(const double multiplier, const Vectorr& rightVector);
	Vectorr operator*(const Vectorr& leftVector, const double multiplier);

	std::ostream& operator<<(std::ostream& outputStream, const Matrixx& outputMatrix);
	std::ostream& operator<<(std::ostream& outputStream, const Roww& outputRow);
	std::ostream& operator<<(std::ostream& outputStream, const Vectorr& outputVector);

	namespace check {
		enum class LengthState {
			NoExcept = 0,
			InvalidHeight = 1,
			InvalidWidth = 2,
			InvalidHeightAndWidth = 3
		};
		struct LengthInfo {
			int height, width;
		};

		enum class IndexState {
			NoExcept = 0,
			RowIndexOutOfRange = 1,
			ColumnIndexOutOfRange = 2,
			BothIndexOutOfRange = 3,
		};
		struct IndexInfo {
			int row, col;
			LengthInfo lengthInfo;
		};

		enum class OperationState {
			NoExcept = 0,
			HeightDoNotMatch = 1,
			WidthDoNotMatch = 2,
			BothLengthDoNotMatch = 3,
			JoinLengthDoNotMatch = 4
		};
		struct OperationInfo {
			char operation;
			LengthInfo lengthInfo1, lengthInfo2;
		};

		const int checkHeight(const int height);
		const int checkWidth(const int width);
		void handleLengthError(const int errorNumber, const LengthInfo lengthInfo);

		const int checkRowIndex(const int row, const int height);
		const int checkColumnIndex(const int col, const int width);
		void handleOutOfRange(const int errorNumber, const IndexInfo indexInfo);

		const int checkHeight(const int height1, const int height2);
		const int checkWidth(const int width1, const int width2);
		const int checkJoinLength(const int width, const int height);
		void handleLogicError(const int errorNumber, const OperationInfo operationInfo);

		double preventNegativeZero(const double value);
	}
}
