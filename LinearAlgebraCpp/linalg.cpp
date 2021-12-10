#include "linalg.h"

namespace linalg {
	Allocator::Allocator(Allocatable& target, const size_t sequence)
		: target(target), sequence(sequence)
	{
	}
	Allocator::~Allocator()
	{
	}
	Allocator& Allocator::operator,(const double value)
	{
		target.allocate(sequence, value);
		sequence++;
		return *this;
	}




	
	Tensor::Tensor()
		: mSize(0)
	{
	}
	Tensor::Tensor(int size)
		: mSize(size)
	{
	}
	/*Tensor::~Tensor()
	{
	}*/

	const size_t Tensor::size() const
	{
		return mSize;
	}
	const double Tensor::convertNegativeZero(const double value) const
	{
		return (value == -0.0) ? 0.0 : value;
	}





	Matrixx::Matrixx()
		: Tensor(), mHeight(0), mWidth(0)
	{
	}
	Matrixx::Matrixx(const int height, const int width)
		: Tensor(height * width)
	{
		init(height, width);
	}
	/*Matrixx::Matrixx(const Matrixx& copyMatrix)
		: Matrixx(static_cast<int>(copyMatrix.mHeight), static_cast<int>(copyMatrix.mWidth))
	{
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				mRows[row][col] = copyMatrix[row][col];
			}
		}
	}*/
	/*Matrixx::Matrixx(Matrixx&& moveMatrix) noexcept
	{
		swap(*this, moveMatrix);
	}*/
	Matrixx::Matrixx(const Roww& copyRow)
		: Matrixx(1, static_cast<int>(copyRow.mSize))
	{
		for (size_t col = 0; col < mHeight; col++) {
			mRows[0][col] = copyRow[col];
		}
	}
	Matrixx::Matrixx(const Vectorr& copyVector)
		: Matrixx(static_cast<int>(copyVector.mSize), 1)
	{
		for (size_t row = 0; row < mHeight; row++) {
			mRows[row][0] = copyVector[row];
		}
	}
	/*Matrixx::~Matrixx()
	{
	}*/
	void Matrixx::init(const int height, const int width)
	{
		int exceptNum = ExceptionHandler::checkValidHeight(height);
		exceptNum += ExceptionHandler::checkValidWidth(width);
		if (exceptNum > static_cast<int>(LengthState::NoExcept)) {
			LengthArgument lengthArg(height, width);
			ExceptionHandler handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}
		
		mHeight = height;
		mWidth = width;

		mRows.clear();
		mRows.resize(mHeight, Roww(width));
	}

	void Matrixx::reduce()
	{
		toEchelonForm();
		toReducedEchelonForm();
	}

	void Matrixx::toEchelonForm()
	{
		size_t beginRow = 0, beginCol = 0;
		while (beginRow < mHeight && beginCol < mWidth) {
			// 1. Find largest absolute value of entries
			Pivot pivot = findPivot(beginRow, beginCol);
			if (pivot.row >= mHeight || pivot.col >= mWidth) {
				// No more pivot
				break;
			}

			// 2. Switch rows to locate pivot into current row
			std::swap(mRows[beginRow], mRows[pivot.row]);
			pivot.row = beginRow;

			// 3. Set zeros under pivot
			replaceRowsUnder(pivot);

			beginRow++;
			beginCol = pivot.col + 1;
		}
	}
	const Matrixx::Pivot Matrixx::findPivot(const size_t beginRow, const size_t beginCol) const
	{
		double maxAbsoluteEntry = 0.0;
		size_t maxAbsoluteRow = mHeight;
		for (size_t col = beginCol; col < mWidth; col++) {
			for (size_t row = beginRow; row < mHeight; row++) {
				if (abs(mRows[row][col]) > maxAbsoluteEntry) {
					maxAbsoluteEntry = abs(mRows[row][col]);
					maxAbsoluteRow = row;
				}
			}
			if (maxAbsoluteRow >= beginRow && maxAbsoluteRow < mHeight) {
				return Pivot{ maxAbsoluteRow, col, mRows[maxAbsoluteRow][col] };
			}
		}
		return Pivot{ mHeight, mWidth, 0.0 }; // Dummy index and value
	}
	const void Matrixx::replaceRowsUnder(const Pivot pivot)
	{
		for (size_t row = pivot.row + 1; row < mHeight; row++) {
			mRows[row] -= (mRows[row][pivot.col] / pivot.entry) * mRows[pivot.row];
		}
	}

	void Matrixx::toReducedEchelonForm()
	{
		if (!isEchelonForm(*this)) {
			EtcArgument etcArg("Cannot reduce non-echelon form matrix.");
			ExceptionHandler handler(ExceptionState::EtcException, static_cast<int>(EtcState::Exception));
			handler.addArgument(etcArg);
			handler.handleException();
		}

		for (int row = static_cast<int>(mHeight) - 1; row >= 0; row--) {
			Pivot pivot = getPivot(row);
			if (pivot.col >= mWidth) {
				// When trying to find pivot in zero rows
				continue;
			}
			// 4. Set zeros over pivot
			replaceRowsOver(pivot);
			// 5. Set pivot as 1
			mRows[row] /= pivot.entry;
		}
	}
	const Matrixx::Pivot Matrixx::getPivot(const size_t row) const
	{
		for (size_t col = 0; col < mWidth; col++) {
			if (mRows[row][col] != 0.0) {
				return Pivot{ row, col, mRows[row][col] };
			}
		}
		return Pivot{ row, mWidth, 0.0}; // Dummy index and value
	}
	const void Matrixx::replaceRowsOver(const Pivot pivot)
	{
		for (size_t row = 0; row < pivot.row; row++) {
			mRows[row] -= (mRows[row][pivot.col] / pivot.entry) * mRows[pivot.row];
		}
	}

	bool isEchelonForm(const Matrixx& matrix)
	{
		Matrixx::Pivot prePivot = matrix.getPivot(0);
		for (size_t row = 1; row < matrix.mHeight; row++) {
			Matrixx::Pivot curPivot = matrix.getPivot(row);
			if (curPivot.col < matrix.mWidth && curPivot.col <= prePivot.col) {
				return false;
			}
			prePivot = curPivot;
		}
		return true;
	}

	Matrixx Matrixx::block(const size_t beginRow, const size_t beginCol,
		const size_t blockHeight, const size_t blockWidth) const
	{
		int exceptNum = ExceptionHandler::checkRowIndex(beginRow, mHeight);
		exceptNum += ExceptionHandler::checkColumnIndex(beginCol, mWidth);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(beginRow, mHeight);
			ColumnIndexArgument colIndexArg(beginCol, mWidth);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}
		exceptNum = ExceptionHandler::checkRowIndex(beginRow + blockHeight - 1, mHeight);
		exceptNum += ExceptionHandler::checkColumnIndex(beginCol + blockWidth - 1, mWidth);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(beginRow + blockHeight - 1, mHeight);
			ColumnIndexArgument colIndexArg(beginCol + blockWidth - 1, mWidth);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		Matrixx blockMatrix(static_cast<int>(blockHeight), static_cast<int>(blockWidth));
		for (size_t row = 0; row < blockMatrix.mHeight; row++) {
			for (size_t col = 0; col < blockMatrix.mWidth; col++) {
				blockMatrix[row][col] = mRows[beginRow + row][beginCol + col];
			}
		}
		return blockMatrix;
	}

	Matrixx Matrixx::inverse()
	{
		if (mHeight != mWidth) {
			EtcArgument etcArg("Cannot get inverse matrix from non-square matrix.");
			ExceptionHandler handler(ExceptionState::EtcException, static_cast<int>(EtcState::Exception));
			handler.addArgument(etcArg);
			handler.handleException();
		}

		const int length = static_cast<int>(mHeight);
		Matrixx identityMatrix = identity(length);
		Matrixx appendedMatrix = *this & identityMatrix;
		appendedMatrix.reduce();

		if (appendedMatrix.block(0, 0, length, length) != identityMatrix) {
			EtcArgument etcArg("The matrix is not reversible.");
			ExceptionHandler handler(ExceptionState::EtcException, static_cast<int>(EtcState::Exception));
			handler.addArgument(etcArg);
			handler.handleException();
		}

		return appendedMatrix.block(0, length, length, length);
	}

	Matrixx Matrixx::transpose(const bool inplace)
	{
		Matrixx transposedMatrix(static_cast<int>(mWidth), static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				transposedMatrix[col][row] = mRows[row][col];
			}
		}
		if (inplace) {
			swap(*this, transposedMatrix);
			return *this;
		}
		else {
			return transposedMatrix;
		}
	}

	Matrixx Matrixx::identity(const int length)
	{
		Matrixx identityMatrix(length, length);
		for (size_t row = 0; row < identityMatrix.mHeight; row++) {
			for (size_t col = 0; col < identityMatrix.mWidth; col++) {
				identityMatrix[row][col] = (row == col) ? 1.0 : 0.0;
			}
		}
		return identityMatrix;
	}

	/*Matrixx Matrixx::zero(const int height, const int width)
	{
		Matrixx zeroMatrix(height, width);
		for (size_t row = 0; row < zeroMatrix.mHeight; row++) {
			for (size_t col = 0; col < zeroMatrix.mWidth; col++) {
				zeroMatrix[row][col] = 0.0;
			}
		}
		return zeroMatrix;
	}*/

	Roww& Matrixx::operator[](const size_t row)
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mHeight);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		return mRows[row];
	}
	const Roww& Matrixx::operator[](const size_t row) const
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mHeight);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		return mRows[row];
	}

	Roww& Matrixx::operator()(const int row)
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mHeight);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		if (row >= 0) {
			return mRows[row];
		}
		else {
			return mRows[static_cast<size_t>(static_cast<int>(mHeight) + row)];
		}
	}
	const Roww& Matrixx::operator()(const int row) const
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mHeight);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		if (row >= 0) {
			return mRows[row];
		}
		else {
			return mRows[static_cast<size_t>(static_cast<int>(mHeight) + row)];
		}
	}

	double& Matrixx::operator()(const int row, const int col)
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mHeight);
		exceptNum += ExceptionHandler::checkColumnIndex(col, mWidth);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight);
			ColumnIndexArgument colIndexArg(col, mWidth);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		if (row >= 0) {
			return mRows[row](col);
		}
		else {
			return mRows[static_cast<size_t>(static_cast<int>(mHeight) + row)](col);
		}
	}
	const double& Matrixx::operator()(const int row, const int col) const
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mHeight);
		exceptNum += ExceptionHandler::checkColumnIndex(col, mWidth);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight);
			ColumnIndexArgument colIndexArg(col, mWidth);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		if (row >= 0) {
			return mRows[row](col);
		}
		else {
			return mRows[static_cast<size_t>(static_cast<int>(mHeight) + row)](col);
		}
	}

	Matrixx& Matrixx::operator=(std::initializer_list<double> values)
	{
		size_t sequence = 0;
		for (double value : values) {
			allocate(sequence, value);
			sequence++;
		}
		return *this;
	}
	Allocator& Matrixx::operator<<(const double value)
	{
		allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Matrixx::allocate(const size_t sequence, const double value)
	{
		if (sequence < mSize) {
			mRows[sequence / mWidth][sequence % mWidth] = convertNegativeZero(value);
		}
	}

	Matrixx Matrixx::operator+() const
	{
		return Matrixx(*this);
	}
	Matrixx Matrixx::operator-() const
	{
		Matrixx negativeMatrix(static_cast<int>(mHeight), static_cast<int>(mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				negativeMatrix[row][col] = convertNegativeZero(-mRows[row][col]);
			}
		}
		return negativeMatrix;
	}

	void swap(Matrixx& leftMatrix, Matrixx& rightMatrix) noexcept
	{
		std::swap(leftMatrix.mRows, rightMatrix.mRows);
		std::swap(leftMatrix.mHeight, rightMatrix.mHeight);
		std::swap(leftMatrix.mWidth, rightMatrix.mWidth);
	}
	/*Matrixx& Matrixx::operator=(const Matrixx& rightMatrix)
	{
		if (this == &rightMatrix) {
			return *this;
		}

		Matrixx copyMatrix(rightMatrix);
		swap(*this, copyMatrix);
		return *this;
	}*/
	/*Matrixx& Matrixx::operator=(Matrixx&& rightMatrix)
	{
		Matrixx moveMatrix(std::move(rightMatrix));
		swap(*this, moveMatrix);
		return *this;
	}*/
	Matrixx& Matrixx::operator+=(const Matrixx& rightMatrix)
	{
		int exceptNum = ExceptionHandler::checkHeight(mHeight, rightMatrix.mHeight);
		exceptNum += ExceptionHandler::checkWidth(mWidth, rightMatrix.mWidth);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightMatrix.mHeight, rightMatrix.mWidth);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Matrixx resultMatrix(*this);
		for (size_t row = 0; row < mHeight; row++) {
			resultMatrix[row] += rightMatrix[row];
		}
		swap(*this, resultMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator-=(const Matrixx& rightMatrix)
	{
		int exceptNum = ExceptionHandler::checkHeight(mHeight, rightMatrix.mHeight);
		exceptNum += ExceptionHandler::checkWidth(mWidth, rightMatrix.mWidth);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightMatrix.mHeight, rightMatrix.mWidth);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Matrixx resultMatrix(*this);
		for (size_t row = 0; row < mHeight; row++) {
			resultMatrix[row] -= rightMatrix[row];
		}
		swap(*this, resultMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator*=(const double multiplier)
	{
		for (size_t row = 0; row < mHeight; row++) {
			mRows[row] *= multiplier;
		}
		return *this;
	}
	Matrixx& Matrixx::operator*=(const Matrixx& rightMatrix)
	{
		int exceptNum = ExceptionHandler::checkJoinLength(mWidth, rightMatrix.mHeight);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightMatrix.mHeight, rightMatrix.mWidth);
			OperationArgument operationArg('*', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Matrixx resultMatrix(static_cast<int>(mHeight), static_cast<int>(rightMatrix.mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < rightMatrix.mWidth; col++) {
				double dotProduct = 0.0;
				for (size_t join = 0; join < mWidth; join++) {
					dotProduct += mRows[row][join] * rightMatrix[join][col];
				}
				resultMatrix[row][col] = convertNegativeZero(dotProduct);
			}
		}
		linalg::swap(*this, resultMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator/=(const double divisor)
	{
		if (convertNegativeZero(divisor) == 0.0) {
			ExceptionHandler handler(ExceptionState::ArithmeticException,
				static_cast<int>(OperationState::DivideByZero));
			handler.handleException();
		}

		Matrixx resultMatrix(*this);
		for (size_t row = 0; row < mHeight; row++) {
			resultMatrix[row] /= divisor;
		}
		swap(*this, resultMatrix);
		return *this;
	}

	Matrixx& Matrixx::operator&=(const Matrixx& rightMatrix)
	{
		int exceptNum = ExceptionHandler::checkHeight(mHeight, rightMatrix.mHeight);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightMatrix.mHeight, rightMatrix.mWidth);
			OperationArgument operationArg('&', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Matrixx appendedMatrix(static_cast<int>(mHeight), static_cast<int>(mWidth + rightMatrix.mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrix[row][col] = mRows[row][col];
			}
		}
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = mWidth; col < mWidth + rightMatrix.mWidth; col++) {
				appendedMatrix[row][col] = rightMatrix[row][col - mWidth];
			}
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator&=(const Vectorr& rightVector)
	{
		int exceptNum = ExceptionHandler::checkHeight(mHeight, rightVector.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightVector.mSize, 1);
			OperationArgument operationArg('&', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Matrixx appendedMatrix(static_cast<int>(mHeight), static_cast<int>(mWidth + 1));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrix[row][col] = mRows[row][col];
			}
		}
		for (size_t row = 0; row < mHeight; row++) {
			appendedMatrix[row][mWidth] = rightVector[row];
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator|=(const Matrixx& lowerMatrix)
	{
		int exceptNum = ExceptionHandler::checkWidth(mWidth, lowerMatrix.mWidth);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument upperLengthArg(mHeight, mWidth);
			LengthArgument lowerLengthArg(lowerMatrix.mHeight, lowerMatrix.mWidth);
			OperationArgument operationArg('|', upperLengthArg, lowerLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Matrixx appendedMatrix(static_cast<int>(mHeight + lowerMatrix.mHeight), static_cast<int>(mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrix[row][col] = mRows[row][col];
			}
		}
		for (size_t row = mHeight; row < mHeight + lowerMatrix.mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrix[row][col] = lowerMatrix[row - mHeight][col];
			}
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator|=(const Roww& lowerRow)
	{
		int exceptNum = ExceptionHandler::checkWidth(mWidth, lowerRow.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument upperLengthArg(mHeight, mWidth);
			LengthArgument lowerLengthArg(1, lowerRow.mSize);
			OperationArgument operationArg('|', upperLengthArg, lowerLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Matrixx appendedMatrix(static_cast<int>(mHeight + 1), static_cast<int>(mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrix[row][col] = mRows[row][col];
			}
		}
		for (size_t col = 0; col < mWidth; col++) {
			appendedMatrix[mHeight][col] = lowerRow[col];
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}

	Matrixx operator+(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		Matrixx resultMatrix(leftMatrix);
		resultMatrix += rightMatrix;
		return resultMatrix;
	}
	
	Matrixx operator-(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		Matrixx resultMatrix(leftMatrix);
		resultMatrix -= rightMatrix;
		return resultMatrix;
	}
	Matrixx operator*(const double multiplier, const Matrixx& rightMatrix)
	{
		Matrixx resultMatrix(rightMatrix);
		resultMatrix *= multiplier;
		return resultMatrix;
	}
	Matrixx operator*(const Matrixx& leftMatrix, const double multiplier)
	{
		return operator*(multiplier, leftMatrix);
	}
	Matrixx operator*(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		Matrixx resultMatrix(leftMatrix);
		resultMatrix *= rightMatrix;
		return resultMatrix;
	}
	Matrixx operator/(const Matrixx& leftMatrix, const double divisor)
	{
		Matrixx resultMatrix(leftMatrix);
		resultMatrix /= divisor;
		return resultMatrix;
	}
	Vectorr operator*(const Matrixx& leftMatrix, const Vectorr& rightVector)
	{
		int exceptNum = ExceptionHandler::checkJoinLength(leftMatrix.mWidth, rightVector.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(leftMatrix.mHeight, leftMatrix.mWidth);
			LengthArgument rightLengthArg(rightVector.mSize, 1);
			OperationArgument operationArg('*', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Vectorr resultVector(static_cast<int>(leftMatrix.mHeight));
		for (size_t row = 0; row < leftMatrix.mHeight; row++) {
			double dotProduct = 0.0;
			for (size_t join = 0; join < leftMatrix.mWidth; join++) {
				dotProduct += leftMatrix[row][join] * rightVector[join];
			}
			resultVector[row] = dotProduct;
		}
		return resultVector;
	}

	Matrixx operator&(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		Matrixx appendedMatrix(leftMatrix);
		appendedMatrix &= rightMatrix;
		return appendedMatrix;
	}
	Matrixx operator&(const Matrixx& leftMatrix, const Vectorr& rightVector)
	{
		Matrixx appendedMatrix(leftMatrix);
		appendedMatrix &= rightVector;
		return appendedMatrix;
	}
	Matrixx operator&(const Vectorr& leftVector, const Matrixx& rightMatrix)
	{
		Matrixx appendedMatrix(leftVector);
		appendedMatrix &= rightMatrix;
		return appendedMatrix;
	}
	Matrixx operator&(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		Matrixx appendedMatrix(leftVector);
		appendedMatrix &= rightVector;
		return appendedMatrix;
	}

	Matrixx operator|(const Matrixx& upperMatrix, const Matrixx& lowerMatrix)
	{
		Matrixx appendedMatrix(upperMatrix);
		appendedMatrix |= lowerMatrix;
		return appendedMatrix;
	}
	Matrixx operator|(const Matrixx& upperMatrix, const Roww& lowerRow)
	{
		Matrixx appendedMatrix(upperMatrix);
		appendedMatrix |= lowerRow;
		return appendedMatrix;
	}
	Matrixx operator|(const Roww& upperRow, const Matrixx& lowerMatrix)
	{
		Matrixx appendedMatrix(upperRow);
		appendedMatrix |= lowerMatrix;
		return appendedMatrix;
	}
	Matrixx operator|(const Roww& upperRow, const Roww& lowerRow)
	{
		Matrixx appendedMatrix(upperRow);
		appendedMatrix |= lowerRow;
		return appendedMatrix;
	}

	bool operator==(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		if (leftMatrix.mHeight != rightMatrix.mHeight ||
			leftMatrix.mWidth != rightMatrix.mWidth) {
			return false;
		}
		for (size_t row = 0; row < leftMatrix.mHeight; row++) {
			if (leftMatrix[row] != rightMatrix[row]) {
				return false;
			}
		}
		return true;
	}
	bool operator!=(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		return !(leftMatrix == rightMatrix);
	}

	std::ostream& operator<<(std::ostream& outputStream, const Matrixx& outputMatrix)
	{
		outputStream << outputMatrix.str();
		return outputStream;
	}
	
	const size_t Matrixx::getHeight() const
	{
		return mHeight;
	}
	const size_t Matrixx::getWidth() const
	{
		return mWidth;
	}

	Roww Matrixx::getRow(const int row) const
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mHeight);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		if (row >= 0) {
			return mRows[row];
		}
		else {
			return mRows[static_cast<size_t>(static_cast<int>(mHeight) + row)];
		}
	}
	Vectorr Matrixx::getColumn(const int col) const
	{
		Vectorr copyVector(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			copyVector[row] = mRows[row](col);
		}
		return copyVector;
	}

	const std::string Matrixx::str() const
	{
		std::string matrixString = "(" + std::to_string(mHeight) + " x " + std::to_string(mWidth) + " Matrix)\n";
		for (size_t row = 0; row < mHeight; row++) {
			matrixString += mRows[row].str();
		}
		return matrixString;
	}






	Roww::Roww()
		: Tensor()
	{
	}
	Roww::Roww(const int size)
		: Tensor(size)
	{
		init(size);
	}
	/*Roww::Roww(const Roww& copyRow)
		: Roww(static_cast<int>(copyRow.mSize))
	{
		for (int col = 0; col < mSize; col++) {
			mEntries[col] = copyRow[col];
		}
	}*/
	/*Roww::Roww(Roww&& moveRow) noexcept
	{
		swap(*this, moveRow);
	}*/
	/*linalg::Roww::~Roww()
	{
	}*/
	void Roww::init(const int size)
	{
		int exceptNum = ExceptionHandler::checkValidWidth(size);
		if (exceptNum > static_cast<int>(LengthState::NoExcept)) {
			LengthArgument lengthArg(1, size);
			ExceptionHandler handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}

		mSize = size;

		mEntries.clear();
		mEntries.resize(mSize, 0.0);
	}

	double& Roww::operator[](const size_t col)
	{
		int exceptNum = ExceptionHandler::checkColumnIndex(col, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			ColumnIndexArgument colIndexArg(col, mSize);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		return mEntries[col];
	}
	const double& Roww::operator[](const size_t col) const
	{
		int exceptNum = ExceptionHandler::checkColumnIndex(col, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			ColumnIndexArgument colIndexArg(col, mSize);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		return mEntries[col];
	}
	double& Roww::operator()(const int col)
	{
		int exceptNum = ExceptionHandler::checkColumnIndex(col, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			ColumnIndexArgument colIndexArg(col, mSize);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		if (col >= 0) {
			return mEntries[col];
		}
		else {
			return mEntries[static_cast<size_t>(static_cast<int>(mSize) + col)];
		}
	}
	const double& Roww::operator()(const int col) const
	{
		int exceptNum = ExceptionHandler::checkColumnIndex(col, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			ColumnIndexArgument colIndexArg(col, mSize);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		if (col >= 0) {
			return mEntries[col];
		}
		else {
			return mEntries[static_cast<size_t>(static_cast<int>(mSize) + col)];
		}
	}

	Roww& Roww::operator=(std::initializer_list<double> values)
	{
		size_t sequence = 0;
		for (double value : values) {
			if (sequence < mSize) {
				allocate(sequence, value);
			}
			sequence++;
		}
		return *this;
	}
	Allocator& Roww::operator<<(const double value)
	{
		allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Roww::allocate(const size_t sequence, const double value)
	{
		if (sequence < mSize) {
			mEntries[sequence] = convertNegativeZero(value);
		}
	}

	Roww Roww::operator+() const
	{
		return Roww(*this);
	}
	Roww Roww::operator-() const
	{
		Roww negativeRow(static_cast<int>(mSize));
		for (size_t col = 0; col < mSize; col++) {
			negativeRow[col] = convertNegativeZero(-mEntries[col]);
		}
		return negativeRow;
	}

	void swap(Roww& leftRow, Roww& rightRow) noexcept
	{
		std::swap(leftRow.mEntries, rightRow.mEntries);
		std::swap(leftRow.mSize, rightRow.mSize);
	}
	/*Roww& Roww::operator=(const Roww& rightRow)
	{
		if (this == &rightRow) {
			return *this;
		}

		Roww copyRow(rightRow);
		swap(*this, copyRow);
		return *this;
	}*/
	/*Roww& Roww::operator=(Roww&& rightRow) noexcept
	{
		Roww moveRow(std::move(rightRow));
		swap(*this, rightRow);
		return *this;
	}*/
	Roww& Roww::operator+=(const Roww& rightRow)
	{
		int exceptNum = ExceptionHandler::checkWidth(mSize, rightRow.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(1, mSize);
			LengthArgument rightLengthArg(1, rightRow.mSize);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Roww resultRow(static_cast<int>(mSize));
		for (size_t col = 0; col < mSize; col++) {
			resultRow[col] = convertNegativeZero(mEntries[col] + rightRow[col]);
		}
		swap(*this, resultRow);
		return *this;
	}
	Roww& Roww::operator-=(const Roww& rightRow)
	{
		int exceptNum = ExceptionHandler::checkWidth(mSize, rightRow.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(1, mSize);
			LengthArgument rightLengthArg(1, rightRow.mSize);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Roww resultRow(static_cast<int>(mSize));
		for (size_t col = 0; col < mSize; col++) {
			resultRow[col] = convertNegativeZero(mEntries[col] - rightRow[col]);
		}
		swap(*this, resultRow);
		return *this;
	}
	Roww& Roww::operator*=(const double multiplier)
	{
		for (size_t col = 0; col < mSize; col++) {
			mEntries[col] = convertNegativeZero(multiplier * mEntries[col]);
		}
		return *this;
	}
	Roww& Roww::operator/=(const double divisor)
	{
		if (convertNegativeZero(divisor) == 0.0) {
			ExceptionHandler handler(ExceptionState::ArithmeticException,
				static_cast<int>(OperationState::DivideByZero));
			handler.handleException();
		}

		Roww resultRow(static_cast<int>(mSize));
		for (size_t col = 0; col < mSize; col++) {
			resultRow[col] = convertNegativeZero(mEntries[col] / divisor);
		}
		swap(*this, resultRow);
		return *this;
	}

	Roww& Roww::operator&=(const Roww& rightRow)
	{
		Roww appendedRow(static_cast<int>(mSize + rightRow.mSize));
		for (size_t col = 0; col < mSize; col++) {
			appendedRow[col] = mEntries[col];
		}
		for (size_t col = mSize; col < mSize + rightRow.mSize; col++) {
			appendedRow[col] = rightRow[col - mSize];
		}
		swap(*this, appendedRow);
		return *this;
	}

	Roww operator+(const Roww& leftRow, const Roww& rightRow)
	{
		Roww resultRow(leftRow);
		resultRow += rightRow;
		return resultRow;
	}
	Roww operator-(const Roww& leftRow, const Roww& rightRow)
	{
		Roww resultRow(leftRow);
		resultRow -= rightRow;
		return resultRow;
	}
	Roww operator*(const double multiplier, const Roww& rightRow)
	{
		Roww resultRow(rightRow);
		resultRow *= multiplier;
		return resultRow;
	}
	Roww operator*(const Roww& leftRow, const double multiplier)
	{
		return operator*(multiplier, leftRow);
	}
	Roww operator/(const Roww& leftRow, const double divisor)
	{
		Roww resultRow(leftRow);
		resultRow /= divisor;
		return resultRow;
	}
	Roww operator&(const Roww& leftRow, const Roww& rightRow)
	{
		Roww appendedRow(leftRow);
		appendedRow &= rightRow;
		return appendedRow;
	}
	bool operator==(const Roww& leftRow, const Roww& rightRow)
	{
		if (leftRow.mSize != rightRow.mSize) {
			return false;
		}
		for (size_t col = 0; col < leftRow.mSize; col++) {
			if (leftRow[col] != rightRow[col]) {
				return false;
			}
		}
		return true;
	}
	bool operator!=(const Roww& leftRow, const Roww& rightRow)
	{
		return !(leftRow == rightRow);
	}
	std::ostream& operator<<(std::ostream& outputStream, const Roww& outputRow)
	{
		outputStream << outputRow.str();
		return outputStream;
	}

	const std::string Roww::str() const
	{
		std::ostringstream parser;
		parser.precision(2);
		std::string rowString = "[\t";
		for (size_t col = 0; col < mSize; col++) {
			if (mEntries[col] >= 0.0) {
				parser << " ";
			}
			parser << mEntries[col];
			rowString += parser.str() + "\t";
			parser.str("");
		}
		rowString += "]\n";
		return rowString;
	}

	
	


	Vectorr::Vectorr()
		: Tensor()
	{
	}
	Vectorr::Vectorr(const int size)
	{
		init(size);
	}
	/*Vectorr::Vectorr(const Vectorr& copyVector)
		: Vectorr(static_cast<int>(copyVector.mSize))
	{
		for (size_t row = 0; row < mSize; row++) {
			mEntries[row] = copyVector[row];
		}
	}*/
	/*Vectorr::Vectorr(Vectorr&& moveVector) noexcept
	{
		swap(*this, moveVector);
	}*/
	/*Vectorr::~Vectorr()
	{
	}*/
	void Vectorr::init(const int size)
	{
		int exceptNum = ExceptionHandler::checkValidHeight(size);
		if (exceptNum > static_cast<int>(LengthState::NoExcept)) {
			LengthArgument lengthArg(size, 1);
			ExceptionHandler handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}

		mSize = size;

		mEntries.clear();
		mEntries.resize(mSize, 0.0);
	}

	double& Vectorr::operator[](const size_t row)
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mSize);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		return mEntries[row];
	}
	const double& Vectorr::operator[](const size_t row) const
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mSize);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		return mEntries[row];
	}
	double& Vectorr::operator()(const int row)
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mSize);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		if (row >= 0) {
			return mEntries[row];
		}
		else {
			return mEntries[static_cast<size_t>(static_cast<int>(mSize) + row)];
		}
	}
	const double& Vectorr::operator()(const int row) const
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mSize);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		if (row >= 0) {
			return mEntries[row];
		}
		else {
			return mEntries[static_cast<size_t>(static_cast<int>(mSize) + row)];
		}
	}

	Vectorr& Vectorr::operator=(std::initializer_list<double> values)
	{
		size_t sequence = 0;
		for (double value : values) {
			allocate(sequence, value);
			sequence++;
		}
		return *this;
	}
	Allocator& Vectorr::operator<<(const double value)
	{
		allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Vectorr::allocate(const size_t sequence, const double value)
	{
		if (sequence < mSize) {
			mEntries[sequence] = convertNegativeZero(value);
		}
	}

	Vectorr Vectorr::operator+() const
	{
		return Vectorr(*this);
	}
	Vectorr Vectorr::operator-() const
	{
		Vectorr negativeVector(static_cast<int>(mSize));
		for (size_t row = 0; row < mSize; row++) {
			negativeVector[row] = convertNegativeZero(-mEntries[row]);
		}
		return negativeVector;
	}

	void swap(Vectorr& leftVector, Vectorr& rightVector) noexcept
	{
		std::swap(leftVector.mEntries, rightVector.mEntries);
		std::swap(leftVector.mSize, rightVector.mSize);
	}
	/*Vectorr& Vectorr::operator=(const Vectorr& rightVector)
	{
		if (this == &rightVector) {
			return *this;
		}

		Vectorr copyVector(rightVector);
		swap(*this, copyVector);
		return *this;
	}*/
	/*Vectorr& Vectorr::operator=(Vectorr&& rightVector) noexcept
	{
		Vectorr moveVector(std::move(rightVector));
		swap(*this, moveVector);
		return *this;
	}*/
	Vectorr& Vectorr::operator+=(const Vectorr& rightVector)
	{
		int exceptNum = ExceptionHandler::checkHeight(mSize, rightVector.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mSize, 1);
			LengthArgument rightLengthArg(rightVector.mSize, 1);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Vectorr resultVector(static_cast<int>(mSize));
		for (size_t row = 0; row < mSize; row++) {
			resultVector[row] = convertNegativeZero(mEntries[row] + rightVector[row]);
		}
		swap(*this, resultVector);
		return *this;
	}
	Vectorr& Vectorr::operator-=(const Vectorr& rightVector)
	{
		int exceptNum = ExceptionHandler::checkHeight(mSize, rightVector.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mSize, 1);
			LengthArgument rightLengthArg(rightVector.mSize, 1);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Vectorr resultVector(static_cast<int>(mSize));
		for (size_t row = 0; row < mSize; row++) {
			resultVector[row] = convertNegativeZero(mEntries[row] - rightVector[row]);
		}
		swap(*this, resultVector);
		return *this;
	}
	Vectorr& Vectorr::operator*=(const double multiplier)
	{
		for (size_t row = 0; row < mSize; row++) {
			mEntries[row] = convertNegativeZero(multiplier * mEntries[row]);
		}
		return *this;
	}
	Vectorr& Vectorr::operator/=(const double divisor)
	{
		if (convertNegativeZero(divisor) == 0.0) {
			ExceptionHandler handler(ExceptionState::ArithmeticException,
				static_cast<int>(OperationState::DivideByZero));
			handler.handleException();
		}

		Vectorr resultVector(static_cast<int>(mSize));
		for (size_t row = 0; row < mSize; row++) {
			resultVector[row] = convertNegativeZero(mEntries[row] / divisor);
		}
		swap(*this, resultVector);
		return *this;
	}

	Vectorr& Vectorr::operator|=(const Vectorr& lowerVector)
	{
		Vectorr appendedVector(static_cast<int>(mSize + lowerVector.mSize));
		for (size_t row = 0; row < mSize; row++) {
			appendedVector[row] = mEntries[row];
		}
		for (size_t row = mSize; row < mSize + lowerVector.mSize; row++) {
			appendedVector[row] = lowerVector[row - mSize];
		}
		swap(*this, appendedVector);
		return *this;
	}

	Vectorr operator+(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		Vectorr resultVector(leftVector);
		resultVector += rightVector;
		return resultVector;
	}
	Vectorr operator-(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		Vectorr resultVector(leftVector);
		resultVector -= rightVector;
		return resultVector;
	}
	Vectorr operator*(const double multiplier, const Vectorr& rightVector)
	{
		Vectorr resultVector(rightVector);
		resultVector *= multiplier;
		return resultVector;
	}
	Vectorr operator*(const Vectorr& leftVector, const double multiplier)
	{
		return operator*(multiplier, leftVector);
	}
	Vectorr operator/(const Vectorr& leftVector, const double divisor)
	{
		Vectorr resultVector(leftVector);
		resultVector /= divisor;
		return resultVector;
	}
	Vectorr operator|(const Vectorr& upperVector, const Vectorr& lowerVector)
	{
		Vectorr appendedVector(upperVector);
		appendedVector |= lowerVector;
		return appendedVector;
	}
	bool operator==(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		if (leftVector.mSize != rightVector.mSize) {
			return false;
		}
		for (size_t row = 0; row < leftVector.mSize; row++) {
			if (leftVector[row] != rightVector[row]) {
				return false;
			}
		}
		return true;
	}
	bool operator!=(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		return !(leftVector == rightVector);
	}
	std::ostream& operator<<(std::ostream& outputStream, const Vectorr& outputVector)
	{
		outputStream << outputVector.str();
		return outputStream;
	}

	const std::string Vectorr::str() const
	{
		std::ostringstream parser;
		parser.precision(2);
		std::string vectorString = "(" + std::to_string(mSize) + " row Vector)\n";
		for (size_t row = 0; row < mSize; row++) {
			if (mEntries[row] >= 0.0) {
				parser << " ";
			}
			parser << mEntries[row];
			vectorString += "[\t" + parser.str() + "\t]\n";
			parser.str("");
		}
		return vectorString;
	}
}