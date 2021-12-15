#include "linalg_exception.h"

namespace linalg {
	ExceptionHandlerr::ExceptionHandlerr(const ExceptionState exceptionState, const int exceptionNumber)
		: mExceptionState(exceptionState), mExceptionNumber(exceptionNumber)
	{
	}
	ExceptionHandlerr::~ExceptionHandlerr()
	{
	}
	void ExceptionHandlerr::setExceptionState(const ExceptionState exceptionState)
	{
		this->mExceptionState = exceptionState;
	}
	void ExceptionHandlerr::setExceptionNumber(const int exceptionNumber)
	{
		this->mExceptionNumber = exceptionNumber;
	}
	void ExceptionHandlerr::addArgument(ExceptionArgument& exceptionArg)
	{
		mExceptionArgs.push_back(&exceptionArg);
	}
	void ExceptionHandlerr::handleException()
	{
		switch (mExceptionState) {
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

	const std::string ExceptionHandlerr::getLengthErrorString() const
	{
		std::string exceptStr;
		switch (mExceptionNumber) {
		case static_cast<int>(LengthState::InvalidHeight):
			exceptStr = "Bad height";
			break;
		case static_cast<int>(LengthState::InvalidWidth):
			exceptStr = "Bad width";
			break;
		case static_cast<int>(LengthState::InvalidHeightAndWidth):
			exceptStr = "Bad height and width";
			break;
		default:
			break;
		}
		if (mExceptionArgs.size() > 0) {
			exceptStr += " : " + mExceptionArgs[0]->str();
		}
		return exceptStr;
	}
	const std::string ExceptionHandlerr::getOutOfRangeString() const
	{
		std::string exceptStr;
		switch (mExceptionNumber) {
		case static_cast<int>(IndexState::RowIndexOutOfRange):
			exceptStr = "Row index out of range.";
			break;
		case static_cast<int>(IndexState::ColumnIndexOutOfRange):
			exceptStr = "Column index out of range.";
			break;
		case static_cast<int>(IndexState::BothIndexOutOfRange):
			exceptStr = "Row and column index out of range.";
			break;
		default:
			break;
		}

		for (size_t argIndex = 0; argIndex < mExceptionArgs.size(); argIndex++) {
			exceptStr += "\n" + mExceptionArgs[argIndex]->str();
		}
		return exceptStr;
	}
	const std::string ExceptionHandlerr::getArithmeticExceptionString() const
	{
		std::string exceptStr;
		switch (mExceptionNumber) {
		case static_cast<int>(OperationState::HeightDoNotMatch):
			exceptStr = "Height do not match";
			break;
		case static_cast<int>(OperationState::WidthDoNotMatch):
			exceptStr = "Width do not match";
			break;
		case static_cast<int>(OperationState::BothLengthDoNotMatch):
			exceptStr = "Height and width do not match";
			break;
		case static_cast<int>(OperationState::JoinLengthDoNotMatch):
			exceptStr = "Cannot multiply";
			break;
		case static_cast<int>(OperationState::DivideByZero):
			return "Divide by zero";
		default:
			break;
		}
		if (mExceptionArgs.size() > 0) {
			exceptStr += " : " + mExceptionArgs[0]->str();
		}
		return exceptStr;
	}
	const std::string ExceptionHandlerr::getEtcExceptionString() const
	{
		if (mExceptionArgs.size() > 0) {
			return mExceptionArgs[0]->str();
		}
		return "Logic Error";
	}

	const int ExceptionHandlerr::checkValidHeight(const size_t height)
	{
		if (height < 1 || height > std::numeric_limits<int>::max()) {
			return static_cast<int>(LengthState::InvalidHeight);
		}
		return static_cast<int>(LengthState::NoExcept);
	}
	const int ExceptionHandlerr::checkValidWidth(const size_t width)
	{
		if (width < 1 || width > std::numeric_limits<int>::max()) {
			return static_cast<int>(LengthState::InvalidWidth);
		}
		return static_cast<int>(LengthState::NoExcept);
	}

	const int ExceptionHandlerr::checkRowIndex(const int row, const size_t height)
	{
		if (row >= static_cast<int>(height) || row < -static_cast<int>(height)) {
			return static_cast<int>(IndexState::RowIndexOutOfRange);
		}
		return static_cast<int>(IndexState::NoExcept);
	}
	const int ExceptionHandlerr::checkRowIndex(const size_t row, const size_t height)
	{
		return ExceptionHandlerr::checkRowIndex(static_cast<int>(row), height);
	}

	const int ExceptionHandlerr::checkColumnIndex(const int col, const size_t width)
	{
		if (col >= static_cast<int>(width) || col < -static_cast<int>(width)) {
			return static_cast<int>(IndexState::ColumnIndexOutOfRange);
		}
		return static_cast<int>(IndexState::NoExcept);
	}
	const int ExceptionHandlerr::checkColumnIndex(const size_t col, const size_t width)
	{
		return ExceptionHandlerr::checkColumnIndex(static_cast<int>(col), width);
	}

	const int ExceptionHandlerr::checkHeight(const size_t height1, const size_t height2)
	{
		if (height1 != height2) {
			return static_cast<int>(OperationState::HeightDoNotMatch);
		}
		return static_cast<int>(OperationState::NoExcept);
	}
	const int ExceptionHandlerr::checkWidth(const size_t width1, const size_t width2)
	{
		if (width1 != width2) {
			return static_cast<int>(OperationState::WidthDoNotMatch);
		}
		return static_cast<int>(OperationState::NoExcept);
	}
	const int ExceptionHandlerr::checkJoinLength(const size_t width, const size_t height)
	{
		if (width != height) {
			return static_cast<int>(OperationState::JoinLengthDoNotMatch);
		}
		return static_cast<int>(OperationState::NoExcept);
	}




	LengthArgument::LengthArgument(const size_t height, const size_t width)
		: mHeight(height), mWidth(width)
	{
	}
	LengthArgument::~LengthArgument()
	{
	}

	const std::string LengthArgument::str() const
	{
		return "Heignt : " + std::to_string(mHeight) + ", Width : " + std::to_string(mWidth);
	}




	RowIndexArgument::RowIndexArgument(const int row, const size_t height, const bool allowNegativeIndex)
		: mRow(row), mHeight(height)
	{
		mBeginIndex = allowNegativeIndex ? -static_cast<int>(height) : 0;
	}
	RowIndexArgument::RowIndexArgument(const size_t row, const size_t height, const bool allowNegativeIndex)
		: RowIndexArgument::RowIndexArgument(static_cast<int>(row), height, allowNegativeIndex)
	{
	}
	RowIndexArgument::~RowIndexArgument()
	{
	}

	const std::string RowIndexArgument::str() const
	{
		return "Row : range(" + std::to_string(mBeginIndex) + " ~ " + std::to_string(mHeight - 1)
			+ "), access(" + std::to_string(mRow) + ")";
	}




	ColumnIndexArgument::ColumnIndexArgument(const int col, const size_t width, const bool allowNegativeIndex)
		: mCol(col), mWidth(width)
	{
		mBeginIndex = allowNegativeIndex ? -static_cast<int>(width) : 0;
	}
	ColumnIndexArgument::ColumnIndexArgument(const size_t col, const size_t width, const bool allowNegativeIndex)
		: ColumnIndexArgument::ColumnIndexArgument(static_cast<int>(col), width, allowNegativeIndex)
	{
	}
	ColumnIndexArgument::~ColumnIndexArgument()
	{
	}

	const std::string ColumnIndexArgument::str() const
	{
		return "Column : range(" + std::to_string(mBeginIndex) + " ~ " + std::to_string(mWidth - 1)
			+ "), access(" + std::to_string(mCol) + ")";
	}




	OperationArgument::OperationArgument(const char operation,
		const LengthArgument& lengthArg1, const LengthArgument& lengthArg2)
		: mOperation(operation), mLengthArg1(lengthArg1), mLengthArg2(lengthArg2)
	{
	}
	OperationArgument::~OperationArgument()
	{
	}

	const std::string OperationArgument::str() const
	{
		return "(" + std::to_string(mLengthArg1.mHeight) + " x " + std::to_string(mLengthArg1.mWidth) + ") "
			+ mOperation
			+ " (" + std::to_string(mLengthArg2.mHeight) + " x " + std::to_string(mLengthArg2.mWidth) + ")";
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