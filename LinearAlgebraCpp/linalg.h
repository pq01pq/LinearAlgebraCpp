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

		Matrixx block(const int beginRow, const int beginCol, const int blockHeight, const int blockWidth) const;

		static Matrixx identity(const int length);
		static Matrixx zero(const int height, const int width);

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
		
		void init(const int width);

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

	

	/*class ExceptionHandler {
	public:
		ExceptionHandler(const ExceptionState exceptionState, const int exceptionNumber);
		~ExceptionHandler();

		void addArgument(const ExceptionArgument& exceptionArg);
		void handleException();

		static const int checkValidHeight(const int height);
		static const int checkValidWidth(const int width);

		static const int checkRowIndex(const int row, const int height);
		static const int checkColumnIndex(const int col, const int width);

		static const int checkHeight(const int height1, const int height2);
		static const int checkWidth(const int width1, const int width2);
		static const int checkJoinLength(const int width, const int height);
	private:
		ExceptionState exceptionState;
		int exceptionNumber;
		std::vector<ExceptionArgument> exceptionArgs;

		const std::string getLengthErrorString() const;
		const std::string getOutOfRangeString() const;
		const std::string getArithmeticExceptionString() const;

		enum class LengthState {
			NoExcept = 0,
			InvalidHeight = 1,
			InvalidWidth = 2,
			InvalidHeightAndWidth = 3
		};
		enum class IndexState {
			NoExcept = 0,
			RowIndexOutOfRange = 1,
			ColumnIndexOutOfRange = 2,
			BothIndexOutOfRange = 3,
		};
		enum class OperationState {
			NoExcept = 0,
			HeightDoNotMatch = 1,
			WidthDoNotMatch = 2,
			BothLengthDoNotMatch = 3,
			JoinLengthDoNotMatch = 4
		};
	};

	enum class ExceptionState {
		NoExcept,
		LengthError,
		OutOfRange,
		ArithmeticException
	};

	class ExceptionArgument {
		friend class ExceptionHandler;
	public:
		virtual const std::string str() const { return ""; };
	};

	class LengthArgument : public ExceptionArgument {
		friend class ExceptionHandler;
		friend class IndexArgument;
		friend class OperationArgument;
	public:
		LengthArgument(const int height, const int width);
		~LengthArgument();

		virtual const std::string str() const override;
	private:
		int height, width;
	};

	class RowIndexArgument : public ExceptionArgument {
		friend class ExceptionHandler;
	public:
		RowIndexArgument(const int row, const int height);
		~RowIndexArgument();

		virtual const std::string str() const override;
	private:
		int row, height;
	};

	class ColumnIndexArgument : public ExceptionArgument {
		friend class ExceptionHandler;
	public:
		ColumnIndexArgument(const int col, const int width);
		~ColumnIndexArgument();

		virtual const std::string str() const override;
	private:
		int col, width;
	};

	class OperationArgument : public ExceptionArgument {
		friend class ExceptionHandler;
	public:
		OperationArgument(char operation, LengthArgument& lengthArg1, LengthArgument& lengthArg2);
		~OperationArgument();

		virtual const std::string str() const override;
	private:
		char operation;
		LengthArgument lengthArg1, lengthArg2;
	};*/
}
