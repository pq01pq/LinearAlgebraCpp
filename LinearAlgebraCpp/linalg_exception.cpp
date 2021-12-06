#include "linalg_exception.h"


namespace linalg {
	ExceptionHandler::ExceptionHandler(const ExceptionState exceptionState, const int exceptionNumber)
		: exceptionState(exceptionState), exceptionNumber(exceptionNumber)
	{
	}
	ExceptionHandler::~ExceptionHandler()
	{
	}
	void ExceptionHandler::setExceptionState(const ExceptionState exceptionState)
	{
		this->exceptionState = exceptionState;
	}
	void ExceptionHandler::setExceptionNumber(const int exceptionNumber)
	{
		this->exceptionNumber = exceptionNumber;
	}
	void ExceptionHandler::addArgument(ExceptionArgument& exceptionArg)
	{
		exceptionArgs.push_back(&exceptionArg);
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
		case ExceptionState::EtcException:
			throw std::logic_error(getEtcExceptionString() + "\n");
		default:
			return;
		}
	}

	const std::string ExceptionHandler::getLengthErrorString() const
	{
		std::string exceptStr;
		switch (exceptionNumber) {
		case (int)LengthState::InvalidHeight:
			exceptStr = "Bad height";
			break;
		case (int)LengthState::InvalidWidth:
			exceptStr = "Bad width";
			break;
		case (int)LengthState::InvalidHeightAndWidth:
			exceptStr = "Bad height and width";
			break;
		default:
			break;
		}
		if (exceptionArgs.size() > 0) {
			exceptStr += " : " + exceptionArgs[0]->str();
		}
		return exceptStr;
	}
	const std::string ExceptionHandler::getOutOfRangeString() const
	{
		std::string exceptStr;
		switch (exceptionNumber) {
		case (int)IndexState::RowIndexOutOfRange:
			exceptStr = "Row index out of range.";
			break;
		case (int)IndexState::ColumnIndexOutOfRange:
			exceptStr = "Column index out of range.";
			break;
		case (int)IndexState::BothIndexOutOfRange:
			exceptStr = "Row and column index out of range.";
			break;
		default:
			break;
		}
		if (exceptionArgs.size() > 0) {
			for (int argIndex = 0; argIndex < exceptionArgs.size(); argIndex++) {
				exceptStr += "\n" + exceptionArgs[argIndex]->str();
			}
		}
		return exceptStr;
	}
	const std::string ExceptionHandler::getArithmeticExceptionString() const
	{
		std::string exceptStr;
		switch (exceptionNumber) {
		case (int)OperationState::HeightDoNotMatch:
			exceptStr = "Height do not match";
			break;
		case (int)OperationState::WidthDoNotMatch:
			exceptStr = "Width do not match";
			break;
		case (int)OperationState::BothLengthDoNotMatch:
			exceptStr = "Height and width do not match";
			break;
		case (int)OperationState::JoinLengthDoNotMatch:
			exceptStr = "Cannot multiply";
			break;
		case (int)OperationState::DivideByZero:
			return "Divide by zero";
		default:
			break;
		}
		if (exceptionArgs.size() > 0) {
			exceptStr += " : " + exceptionArgs[0]->str();
		}
		return exceptStr;
	}
	const std::string ExceptionHandler::getEtcExceptionString() const
	{
		if (exceptionArgs.size() > 0) {
			return exceptionArgs[0]->str();
		}
		return "Logic Error";
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
		if (row >= height || row < -height) {
			return (int)IndexState::RowIndexOutOfRange;
		}
		return (int)IndexState::NoExcept;
	}
	const int ExceptionHandler::checkColumnIndex(const int col, const int width)
	{
		if (col >= width || col < -width) {
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
		return "Row : range(" + std::to_string(-height) + " ~ " + std::to_string(height - 1)
			+ "), access(" + std::to_string(row) + ")";
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
		return "Column : range(" + std::to_string(-width) + " ~ " + std::to_string(width - 1)
			+ "), access(" + std::to_string(col) + ")";
	}




	OperationArgument::OperationArgument(const char operation,
		const LengthArgument& lengthArg1, const LengthArgument& lengthArg2)
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




	EtcArgument::EtcArgument(const std::string& what)
		: what(what)
	{
	}
	EtcArgument::~EtcArgument()
	{
	}
	const std::string EtcArgument::str() const
	{
		return what;
	}
}