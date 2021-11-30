#pragma once
#include <string>
#include <vector>

namespace linalg {

	class ExceptionHandler;
	class ExceptionArgument;
	class LengthArgument;
	class IndexArgument;
	class OperationArgument;
	enum class ExceptionState;

	class ExceptionHandler {
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
	};
}