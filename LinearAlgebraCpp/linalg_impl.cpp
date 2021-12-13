#include "linalg_impl.h"

namespace linalg {
	Matrixx::Impl::Impl(const int height, const int width)
	{
		init(height, width);
	}
	/*Matrixx::Matrixx(Matrixx&& moveMatrix) noexcept
	{
		swap(*this, moveMatrix);
	}*/
	Matrixx::Impl::Impl(const Roww& copyRow)
		: Impl(1, static_cast<int>(copyRow.width()))
	{
		for (size_t col = 0; col < mHeight; col++) {
			mRows[0][col] = copyRow[col];
		}
	}
	Matrixx::Impl::Impl(const Vectorr& copyVector)
		: Impl(static_cast<int>(copyVector.height()), 1)
	{
		for (size_t row = 0; row < mHeight; row++) {
			mRows[row][0] = copyVector[row];
		}
	}
	void Matrixx::Impl::init(const int height, const int width)
	{
		int exceptNum = ExceptionHandlerr::checkValidHeight(height);
		exceptNum += ExceptionHandlerr::checkValidWidth(width);
		if (exceptNum > static_cast<int>(LengthState::NoExcept)) {
			LengthArgument lengthArg(height, width);
			ExceptionHandlerr handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}

		mHeight = height;
		mWidth = width;

		mRows.clear();
		mRows.resize(mHeight, Roww(width));
	}

	void Matrixx::Impl::reduce()
	{
		toEchelonForm();
		toReducedEchelonForm();
	}

	void Matrixx::Impl::toEchelonForm()
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
	const Matrixx::Impl::Pivot Matrixx::Impl::findPivot(const size_t beginRow, const size_t beginCol) const
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
	const void Matrixx::Impl::replaceRowsUnder(const Pivot pivot)
	{
		for (size_t row = pivot.row + 1; row < mHeight; row++) {
			mRows[row] -= (mRows[row][pivot.col] / pivot.entry) * mRows[pivot.row];
		}
	}

	void Matrixx::Impl::toReducedEchelonForm()
	{
		if (!this->isEchelonForm()) {
			EtcArgument etcArg("Cannot reduce non-echelon form matrix.");
			ExceptionHandlerr handler(ExceptionState::EtcException, static_cast<int>(EtcState::Exception));
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
	const Matrixx::Impl::Pivot Matrixx::Impl::getPivot(const size_t row) const
	{
		for (size_t col = 0; col < mWidth; col++) {
			if (mRows[row][col] != 0.0) {
				return Pivot{ row, col, mRows[row][col] };
			}
		}
		return Pivot{ row, mWidth, 0.0 }; // Dummy index and value
	}
	const void Matrixx::Impl::replaceRowsOver(const Pivot pivot)
	{
		for (size_t row = 0; row < pivot.row; row++) {
			mRows[row] -= (mRows[row][pivot.col] / pivot.entry) * mRows[pivot.row];
		}
	}

	bool Matrixx::Impl::isEchelonForm()
	{
		Pivot prePivot = getPivot(0);
		for (size_t row = 1; row < mHeight; row++) {
			Pivot curPivot = getPivot(row);
			if (curPivot.col < mWidth && curPivot.col <= prePivot.col) {
				return false;
			}
			prePivot = curPivot;
		}
		return true;
	}

	Matrixx Matrixx::Impl::block(const size_t beginRow, const size_t beginCol,
		const size_t blockHeight, const size_t blockWidth) const
	{
		int exceptNum = ExceptionHandlerr::checkRowIndex(beginRow, mHeight);
		exceptNum += ExceptionHandlerr::checkColumnIndex(beginCol, mWidth);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(beginRow, mHeight);
			ColumnIndexArgument colIndexArg(beginCol, mWidth);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}
		exceptNum = ExceptionHandlerr::checkRowIndex(beginRow + blockHeight - 1, mHeight);
		exceptNum += ExceptionHandlerr::checkColumnIndex(beginCol + blockWidth - 1, mWidth);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(beginRow + blockHeight - 1, mHeight);
			ColumnIndexArgument colIndexArg(beginCol + blockWidth - 1, mWidth);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		Matrixx blockMatrix(static_cast<int>(blockHeight), static_cast<int>(blockWidth));
		for (size_t row = 0; row < blockMatrix.height(); row++) {
			for (size_t col = 0; col < blockMatrix.width(); col++) {
				blockMatrix[row][col] = mRows[beginRow + row][beginCol + col];
			}
		}
		return blockMatrix;
	}

	Matrixx Matrixx::Impl::inverse()
	{
		if (mHeight != mWidth) {
			EtcArgument etcArg("Cannot get inverse matrix from non-square matrix.");
			ExceptionHandlerr handler(ExceptionState::EtcException, static_cast<int>(EtcState::Exception));
			handler.addArgument(etcArg);
			handler.handleException();
		}

		const int length = static_cast<int>(mHeight);
		Matrixx identityMatrix = identity(length);
		Matrixx appendedMatrix = matrix(*this) & identityMatrix;
		appendedMatrix.reduce();

		if (appendedMatrix.block(0, 0, length, length) != identityMatrix) {
			EtcArgument etcArg("The matrix is not reversible.");
			ExceptionHandlerr handler(ExceptionState::EtcException, static_cast<int>(EtcState::Exception));
			handler.addArgument(etcArg);
			handler.handleException();
		}

		return appendedMatrix.block(0, length, length, length);
	}

	Matrixx Matrixx::Impl::transpose() const
	{
		Matrixx transposedMatrix(static_cast<int>(mWidth), static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				transposedMatrix[col][row] = mRows[row][col];
			}
		}
		return transposedMatrix;
	}

	Matrixx Matrixx::Impl::identity(const int length)
	{
		Matrixx identityMatrix(length, length);
		for (size_t row = 0; row < identityMatrix.height(); row++) {
			for (size_t col = 0; col < identityMatrix.width(); col++) {
				identityMatrix[row][col] = (row == col) ? 1.0 : 0.0;
			}
		}
		return identityMatrix;
	}

	Matrixx Matrixx::Impl::matrix(const Impl& matrixImpl)
	{
		Matrixx matrix(static_cast<int>(matrixImpl.mHeight), static_cast<int>(matrixImpl.mWidth));
		for (size_t row = 0; row < matrix.height(); row++) {
			for (size_t col = 0; col < matrix.width(); col++) {
				matrix[row][col] = matrixImpl.get(static_cast<int>(row), static_cast<int>(col));
			}
		}
		return matrix;
	}

	void Matrixx::Impl::allocate(const size_t sequence, const double value)
	{
		if (sequence < mHeight * mWidth) {
			mRows[sequence / mWidth][sequence % mWidth] = convertNegativeZero(value);
		}
	}
	void Matrixx::Impl::allocate(const std::initializer_list<double> values)
	{
		size_t sequence = 0;
		for (double value : values) {
			allocate(sequence, value);
			sequence++;
		}
	}

	Matrixx Matrixx::Impl::negative() const
	{
		Matrixx negativeMatrix(static_cast<int>(mHeight), static_cast<int>(mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				negativeMatrix[row][col] = convertNegativeZero(-mRows[row][col]);
			}
		}
		return negativeMatrix;
	}

	

	void Matrixx::Impl::swap(Impl& rightMatrixImpl) noexcept
	{
		std::swap(mRows, rightMatrixImpl.mRows);
		std::swap(mHeight, rightMatrixImpl.mHeight);
		std::swap(mWidth, rightMatrixImpl.mWidth);
	}

	void Matrixx::Impl::innerAdd(const Matrixx& rightMatrix)
	{
		Impl resultMatrixImpl(*this);
		for (size_t row = 0; row < mHeight; row++) {
			resultMatrixImpl.row(row) += rightMatrix[row];
		}
		swap(resultMatrixImpl);
	}
	void Matrixx::Impl::innerSub(const Matrixx& rightMatrix)
	{
		Impl resultMatrixImpl(*this);
		for (size_t row = 0; row < mHeight; row++) {
			resultMatrixImpl.row(row) -= rightMatrix[row];
		}
		swap(resultMatrixImpl);
	}
	void Matrixx::Impl::innerMul(const double multiplier)
	{
		for (size_t row = 0; row < mHeight; row++) {
			mRows[row] *= multiplier;
		}
	}
	void Matrixx::Impl::innerMatMul(const Matrixx& rightMatrix)
	{
		Impl resultMatrixImpl(static_cast<int>(mHeight), static_cast<int>(rightMatrix.width()));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < rightMatrix.width(); col++) {
				double dotProduct = 0.0;
				for (size_t join = 0; join < mWidth; join++) {
					dotProduct += mRows[row][join] * rightMatrix[join][col];
				}
				resultMatrixImpl.row(row)[col] = convertNegativeZero(dotProduct);
			}
		}
		swap(resultMatrixImpl);
	}
	void Matrixx::Impl::innerDiv(const double divisor)
	{
		Impl resultMatrixImpl(*this);
		for (size_t row = 0; row < mHeight; row++) {
			resultMatrixImpl.row(row) /= divisor;
		}
		swap(resultMatrixImpl);
	}
	void Matrixx::Impl::innerHorizontalAppend(const Matrixx& rightMatrix)
	{
		Impl appendedMatrixImpl(static_cast<int>(mHeight), static_cast<int>(mWidth + rightMatrix.width()));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrixImpl.row(row)[col] = mRows[row][col];
			}
		}
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = mWidth; col < mWidth + rightMatrix.width(); col++) {
				appendedMatrixImpl.row(row)[col] = rightMatrix[row][col - mWidth];
			}
		}
		swap(appendedMatrixImpl);
	}
	void Matrixx::Impl::innerHorizontalAppend(const Vectorr& rightVector)
	{
		innerHorizontalAppend(Matrixx(rightVector));
	}
	void Matrixx::Impl::innerVerticalAppend(const Matrixx& lowerMatrix)
	{
		Impl appendedMatrixImpl(static_cast<int>(mHeight + lowerMatrix.height()), static_cast<int>(mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrixImpl.row(row)[col] = mRows[row][col];
			}
		}
		for (size_t row = mHeight; row < mHeight + lowerMatrix.height(); row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrixImpl.row(row)[col] = lowerMatrix[row - mHeight][col];
			}
		}
		swap(appendedMatrixImpl);
	}
	void Matrixx::Impl::innerVerticalAppend(const Roww& lowerRow)
	{
		innerVerticalAppend(Matrixx(lowerRow));
	}

	Vectorr Matrixx::Impl::vectorEquation(const Vectorr& rightVector)
	{
		Vectorr resultVector(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			double dotProduct = 0.0;
			for (size_t join = 0; join < mWidth; join++) {
				dotProduct += mRows[row][join] * rightVector[join];
			}
			resultVector[row] = dotProduct;
		}
		return resultVector;
	}

	bool Matrixx::Impl::equals(const Matrixx& rightMatrix) const
	{
		if (mHeight != rightMatrix.height() ||
			mWidth != rightMatrix.width()) {
			return false;
		}
		for (size_t row = 0; row < mHeight; row++) {
			if (mRows[row] != rightMatrix[row]) {
				return false;
			}
		}
		return true;
	}

	const Roww& Matrixx::Impl::row(const size_t row) const
	{
		return mRows[row];
	}
	Roww& Matrixx::Impl::row(const size_t row)
	{
		return const_cast<Roww&>(static_cast<const Impl&>(*this).row(row));
	}
	const Roww& Matrixx::Impl::row(const int row) const
	{
		if (row >= 0) {
			return mRows[row];
		}
		else {
			return mRows[static_cast<size_t>(static_cast<int>(mHeight) + row)];
		}
	}
	Roww& Matrixx::Impl::row(const int row)
	{
		return const_cast<Roww&>(static_cast<const Impl&>(*this).row(row));
	}
	const double& Matrixx::Impl::get(const int row, const int col) const
	{
		if (row >= 0) {
			return mRows[row](col);
		}
		else {
			return mRows[static_cast<size_t>(static_cast<int>(mHeight) + row)](col);
		}
	}
	double& Matrixx::Impl::get(const int row, const int col)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this).get(row, col));
	}

	const Roww Matrixx::Impl::getRow(const int row) const
	{
		if (row >= 0) {
			return mRows[row];
		}
		else {
			return mRows[static_cast<size_t>(static_cast<int>(mHeight) + row)];
		}
	}
	const Vectorr Matrixx::Impl::getColumn(const int col) const
	{
		Vectorr copyVector(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			copyVector[row] = mRows[row](col);
		}
		return copyVector;
	}

	const size_t Matrixx::Impl::height() const
	{
		return mHeight;
	}
	const size_t Matrixx::Impl::width() const
	{
		return mWidth;
	}

	const std::string Matrixx::Impl::str() const
	{
		std::string matrixString = "(" + std::to_string(mHeight) + " x " + std::to_string(mWidth) + " Matrix)\n";
		for (size_t row = 0; row < mHeight; row++) {
			matrixString += mRows[row].str();
		}
		return matrixString;
	}





	Roww::Impl::Impl(const int size)
	{
		init(size);
	}
	void Roww::Impl::init(const int size)
	{
		int exceptNum = ExceptionHandlerr::checkValidWidth(size);
		if (exceptNum > static_cast<int>(LengthState::NoExcept)) {
			LengthArgument lengthArg(1, size);
			ExceptionHandlerr handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}

		mWidth = size;

		mEntries.clear();
		mEntries.resize(mWidth, 0.0);
	}

	void Roww::Impl::allocate(const size_t sequence, const double value)
	{
		if (sequence < mWidth) {
			mEntries[sequence] = convertNegativeZero(value);
		}
	}
	void Roww::Impl::allocate(const std::initializer_list<double> values)
	{
		size_t sequence = 0;
		for (double value : values) {
			allocate(sequence, value);
			sequence++;
		}
	}

	Roww Roww::Impl::negative() const
	{
		Roww negativeRow(static_cast<int>(mWidth));
		for (size_t col = 0; col < mWidth; col++) {
			negativeRow[col] = convertNegativeZero(-mEntries[col]);
		}
		return negativeRow;
	}

	void Roww::Impl::swap(Impl& rightRowImpl) noexcept
	{
		std::swap(mEntries, rightRowImpl.mEntries);
		std::swap(mWidth, rightRowImpl.mWidth);
	}

	void Roww::Impl::innerAdd(const Roww& rightRow)
	{
		Impl resultRowImpl(static_cast<int>(mWidth));
		for (size_t col = 0; col < mWidth; col++) {
			resultRowImpl.get(col) = convertNegativeZero(mEntries[col] + rightRow[col]);
		}
		swap(resultRowImpl);
	}
	void Roww::Impl::innerSub(const Roww& rightRow)
	{
		Impl resultRowImpl(static_cast<int>(mWidth));
		for (size_t col = 0; col < mWidth; col++) {
			resultRowImpl.get(col) = convertNegativeZero(mEntries[col] - rightRow[col]);
		}
		swap(resultRowImpl);
	}
	void Roww::Impl::innerMul(const double multiplier)
	{
		for (size_t col = 0; col < mWidth; col++) {
			mEntries[col] = convertNegativeZero(multiplier * mEntries[col]);
		}
	}
	void Roww::Impl::innerDiv(const double divisor)
	{
		Impl resultRowImpl(static_cast<int>(mWidth));
		for (size_t col = 0; col < mWidth; col++) {
			resultRowImpl.get(col) = convertNegativeZero(mEntries[col] / divisor);
		}
		swap(resultRowImpl);
	}
	void Roww::Impl::innerHorizontalAppend(const Roww& rightRow)
	{
		Impl appendedRowImpl(static_cast<int>(mWidth + rightRow.width()));
		for (size_t col = 0; col < mWidth; col++) {
			appendedRowImpl.get(col) = mEntries[col];
		}
		for (size_t col = mWidth; col < mWidth + rightRow.width(); col++) {
			appendedRowImpl.get(col) = rightRow[col - mWidth];
		}
		swap(appendedRowImpl);
	}

	bool Roww::Impl::equals(const Roww& rightRow) const
	{
		if (mWidth != rightRow.width()) {
			return false;
		}
		for (size_t col = 0; col < mWidth; col++) {
			if (mEntries[col] != rightRow[col]) {
				return false;
			}
		}
		return true;
	}

	const double& Roww::Impl::get(const size_t col) const
	{
		return mEntries[col];
	}
	double& Roww::Impl::get(const size_t col)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this).get(col));
	}
	const double& Roww::Impl::get(const int col) const
	{
		if (col >= 0) {
			return mEntries[col];
		}
		else {
			return mEntries[static_cast<size_t>(static_cast<int>(mWidth) + col)];
		}
	}
	double& Roww::Impl::get(const int col)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this).get(col));
	}

	const size_t Roww::Impl::width() const
	{
		return mWidth;
	}

	const std::string Roww::Impl::str() const
	{
		std::ostringstream parser;
		parser.precision(2);
		std::string rowString = "[\t";
		for (size_t col = 0; col < mWidth; col++) {
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





	Vectorr::Impl::Impl(const int size)
	{
		init(size);
	}
	void Vectorr::Impl::init(const int size)
	{
		int exceptNum = ExceptionHandlerr::checkValidHeight(size);
		if (exceptNum > static_cast<int>(LengthState::NoExcept)) {
			LengthArgument lengthArg(size, 1);
			ExceptionHandlerr handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}

		mHeight = size;

		mEntries.clear();
		mEntries.resize(mHeight, 0.0);
	}

	void Vectorr::Impl::allocate(const size_t sequence, const double value)
	{
		if (sequence < mHeight) {
			mEntries[sequence] = convertNegativeZero(value);
		}
	}
	void Vectorr::Impl::allocate(const std::initializer_list<double> values)
	{
		size_t sequence = 0;
		for (double value : values) {
			allocate(sequence, value);
			sequence++;
		}
	}

	Vectorr Vectorr::Impl::negative() const
	{
		Vectorr negativeVector(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			negativeVector[row] = convertNegativeZero(-mEntries[row]);
		}
		return negativeVector;
	}

	void Vectorr::Impl::swap(Impl& rightVectorImpl) noexcept
	{
		std::swap(mEntries, rightVectorImpl.mEntries);
		std::swap(mHeight, rightVectorImpl.mHeight);
	}

	void Vectorr::Impl::innerAdd(const Vectorr& rightVector)
	{
		Impl resultVectorImpl(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			resultVectorImpl.get(row) = convertNegativeZero(mEntries[row] + rightVector[row]);
		}
		swap(resultVectorImpl);
	}

	void Vectorr::Impl::innerSub(const Vectorr& rightVector)
	{
		Impl resultVectorImpl(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			resultVectorImpl.get(row) = convertNegativeZero(mEntries[row] - rightVector[row]);
		}
		swap(resultVectorImpl);
	}
	void Vectorr::Impl::innerMul(const double multiplier)
	{
		for (size_t row = 0; row < mHeight; row++) {
			mEntries[row] = convertNegativeZero(multiplier * mEntries[row]);
		}
	}
	void Vectorr::Impl::innerDiv(const double divisor)
	{
		Impl resultVectorImpl(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			resultVectorImpl.get(row) = convertNegativeZero(mEntries[row] / divisor);
		}
		swap(resultVectorImpl);
	}
	void Vectorr::Impl::innerVerticalAppend(const Vectorr& lowerVector)
	{
		Impl appendedVectorImpl(static_cast<int>(mHeight + lowerVector.height()));
		for (size_t row = 0; row < mHeight; row++) {
			appendedVectorImpl.get(row) = mEntries[row];
		}
		for (size_t row = mHeight; row < mHeight + lowerVector.height(); row++) {
			appendedVectorImpl.get(row) = lowerVector[row - mHeight];
		}
		swap(appendedVectorImpl);
	}

	bool Vectorr::Impl::equals(const Vectorr& rightVector) const
	{
		if (mHeight != rightVector.height()) {
			return false;
		}
		for (size_t row = 0; row < mHeight; row++) {
			if (mEntries[row] != rightVector[row]) {
				return false;
			}
		}
		return true;
	}

	const double& Vectorr::Impl::get(const size_t row) const
	{
		return mEntries[row];
	}
	double& Vectorr::Impl::get(const size_t row)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this).get(row));
	}
	const double& Vectorr::Impl::get(const int row) const
	{
		if (row >= 0) {
			return mEntries[row];
		}
		else {
			return mEntries[static_cast<size_t>(static_cast<int>(mHeight) + row)];
		}
	}
	double& Vectorr::Impl::get(const int row)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this).get(row));
	}

	const size_t Vectorr::Impl::height() const
	{
		return mHeight;
	}

	const std::string Vectorr::Impl::str() const
	{
		std::ostringstream parser;
		parser.precision(2);
		std::string vectorString = "(" + std::to_string(mHeight) + " row Vector)\n";
		for (size_t row = 0; row < mHeight; row++) {
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