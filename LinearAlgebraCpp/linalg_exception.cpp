#include "linalg_exception.h"

namespace linalg {
	ExceptionHandler::ExceptionHandler(const ExceptionState exceptionState, const int exceptionNumber)
		: exceptionState(exceptionState), exceptionNumber(exceptionNumber)
	{
	}
	ExceptionHandler::~ExceptionHandler()
	{
	}
	void ExceptionHandler::addArgument(const ExceptionArgument& exceptionArg)
	{
		exceptionArgs.push_back(exceptionArg);
	}
	void ExceptionHandler::handleException()
	{
		switch (exceptionState) {
		case ExceptionState::LengthError:
			throw std::length_error(getLengthErrorString() + "\n");
		case ExceptionState::OutOfRange:
			throw std::out_of_range(getOutOfRangeString() + "\n");
		case ExceptionState::ArithmeticException:
			throw std::logic_error(getArithmeticExceptionString() + "\n");
		default:
			return;
		}
	}

	const std::string ExceptionHandler::getLengthErrorString() const
	{
		std::string exceptionString;
		switch (exceptionNumber) {
		case (int)LengthState::InvalidHeight:
			exceptionString = "Bad height";
		case (int)LengthState::InvalidWidth:
			exceptionString = "Bad width";
		case (int)LengthState::InvalidHeightAndWidth:
			exceptionString = "Bad height and width";
		default:
			break;
		}
		if (exceptionArgs.size() > 0) {
			exceptionString += " : " + exceptionArgs[0].str();
		}
		return exceptionString;
	}
	const std::string ExceptionHandler::getOutOfRangeString() const
	{
		std::string exceptionString;
		switch (exceptionNumber) {
		case (int)IndexState::RowIndexOutOfRange:
			exceptionString = "Row index out of range";
		case (int)IndexState::ColumnIndexOutOfRange:
			exceptionString = "Column index out of range";
		case (int)IndexState::BothIndexOutOfRange:
			exceptionString = "Row and column index out of range";
		default:
			break;
		}
		if (exceptionArgs.size() > 0) {
			for (int argIndex = 0; argIndex < exceptionArgs.size(); argIndex++) {
				exceptionString += "\n" + exceptionArgs[argIndex].str();
			}
		}
		return exceptionString;
	}
	const std::string ExceptionHandler::getArithmeticExceptionString() const
	{
		std::string exceptionString;
		switch (exceptionNumber) {
		case (int)OperationState::HeightDoNotMatch:
			exceptionString = "Height do not match";
		case (int)OperationState::WidthDoNotMatch:
			exceptionString = "Width do not match";
		case (int)OperationState::BothLengthDoNotMatch:
			exceptionString = "Height and width do not match";
		case (int)OperationState::JoinLengthDoNotMatch:
			exceptionString = "Cannot multiply";
		default:
			break;
		}
		if (exceptionArgs.size() > 0) {
			exceptionString += " : " + exceptionArgs[0].str();
		}
		return exceptionString;
	}

	const int ExceptionHandler::checkValidHeight(const int height)
	{
		if (height < 1) {
			return (int)LengthState::InvalidHeight;
		}
		return (int)LengthState::NoExcept;
	}
	const int ExceptionHandler::checkValidWidth(const int width)
	{
		if (width < 1) {
			return (int)LengthState::InvalidWidth;
		}
		return (int)LengthState::NoExcept;
	}

	const int ExceptionHandler::checkRowIndex(const int row, const int height)
	{
		if (row >= height || row < 0) {
			return (int)IndexState::RowIndexOutOfRange;
		}
		return (int)IndexState::NoExcept;
	}
	const int ExceptionHandler::checkColumnIndex(const int col, const int width)
	{
		if (col >= width || col < 0) {
			return (int)IndexState::ColumnIndexOutOfRange;
		}
		return (int)IndexState::NoExcept;
	}

	const int ExceptionHandler::checkHeight(const int height1, const int height2)
	{
		if (height1 != height2) {
			return (int)OperationState::HeightDoNotMatch;
		}
		return (int)OperationState::NoExcept;
	}
	const int ExceptionHandler::checkWidth(const int width1, const int width2)
	{
		if (width1 != width2) {
			return (int)OperationState::WidthDoNotMatch;
		}
		return (int)OperationState::NoExcept;
	}
	const int ExceptionHandler::checkJoinLength(const int width, const int height)
	{
		if (width != height) {
			return (int)OperationState::JoinLengthDoNotMatch;
		}
		return (int)OperationState::NoExcept;
	}





	LengthArgument::LengthArgument(const int height, const int width)
		: height(height), width(width)
	{
	}
	LengthArgument::~LengthArgument()
	{
	}

	const std::string LengthArgument::str() const
	{
		return "Heignt : " + std::to_string(height) + ", Width : " + std::to_string(width);
	}




	RowIndexArgument::RowIndexArgument(const int row, const int height)
		: row(row), height(height)
	{
	}
	RowIndexArgument::~RowIndexArgument()
	{
	}

	const std::string RowIndexArgument::str() const
	{
		return "Row : range(0-" + std::to_string(row - 1) + "), access(" + std::to_string(height) + ")";
	}



	ColumnIndexArgument::ColumnIndexArgument(const int col, const int width)
		: col(col), width(width)
	{
	}
	ColumnIndexArgument::~ColumnIndexArgument()
	{
	}

	const std::string ColumnIndexArgument::str() const
	{
		return "Column : range(0-" + std::to_string(col - 1) + "), access(" + std::to_string(width) + ")";
	}



	OperationArgument::OperationArgument(char operation, LengthArgument& lengthArg1, LengthArgument& lengthArg2)
		: operation(operation), lengthArg1(lengthArg1), lengthArg2(lengthArg2)
	{
	}
	OperationArgument::~OperationArgument()
	{
	}

	const std::string OperationArgument::str() const
	{
		return "(" + std::to_string(lengthArg1.height) + " x " + std::to_string(lengthArg1.width) + ") "
			+ operation
			+ " (" + std::to_string(lengthArg2.height) + " x " + std::to_string(lengthArg2.width) + ")";
	}
}