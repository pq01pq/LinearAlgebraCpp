#include "linalg_impl.h"

namespace linalg {

	const double Tensorr::Impl::convertNegativeZero(const double value)
	{
		return (value == -0.0) ? 0.0 : value;
	}

	Matrixx::Impl::Impl(const int height, const int width)
	{
		init(height, width);
	}
	Matrixx::Impl::Impl(const Roww::Impl& copyRowImpl)
		: Impl(1, static_cast<int>(copyRowImpl.width()))
	{
		for (size_t col = 0; col < mHeight; col++) {
			mRows[0][col] = copyRowImpl[col];
		}
	}
	Matrixx::Impl::Impl(const Vectorr::Impl& copyVectorImpl)
		: Impl(static_cast<int>(copyVectorImpl.height()), 1)
	{
		for (size_t row = 0; row < mHeight; row++) {
			mRows[row][0] = copyVectorImpl[row];
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
		mRows.resize(mHeight, Roww::Impl(width));
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
			mRows[row] -= mRows[pivot.row] * (mRows[row][pivot.col] / pivot.entry);
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
			mRows[row] -= mRows[pivot.row] * (mRows[row][pivot.col] / pivot.entry);
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

	Matrixx::Impl Matrixx::Impl::block(const size_t beginRow, const size_t beginCol,
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

		Impl blockMatrixImpl(static_cast<int>(blockHeight), static_cast<int>(blockWidth));
		for (size_t row = 0; row < blockMatrixImpl.height(); row++) {
			for (size_t col = 0; col < blockMatrixImpl.width(); col++) {
				blockMatrixImpl[row][col] = mRows[beginRow + row][beginCol + col];
			}
		}
		return blockMatrixImpl;
	}

	Matrixx::Impl Matrixx::Impl::inverse()
	{
		if (mHeight != mWidth) {
			EtcArgument etcArg("Cannot get inverse matrix from non-square matrix.");
			ExceptionHandlerr handler(ExceptionState::EtcException, static_cast<int>(EtcState::Exception));
			handler.addArgument(etcArg);
			handler.handleException();
		}

		const int length = static_cast<int>(mHeight);
		Impl identityMatrixImpl = identity(length);
		Impl appendedMatrixImpl = *this & identityMatrixImpl;
		appendedMatrixImpl.reduce();

		if (appendedMatrixImpl.block(0, 0, length, length) != identityMatrixImpl) {
			EtcArgument etcArg("The matrix is not reversible.");
			ExceptionHandlerr handler(ExceptionState::EtcException, static_cast<int>(EtcState::Exception));
			handler.addArgument(etcArg);
			handler.handleException();
		}

		return appendedMatrixImpl.block(0, length, length, length);
	}

	Matrixx::Impl Matrixx::Impl::transpose() const
	{
		Impl transposedMatrixImpl(static_cast<int>(mWidth), static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				transposedMatrixImpl[col][row] = mRows[row][col];
			}
		}
		return transposedMatrixImpl;
	}

	Matrixx::Impl Matrixx::Impl::identity(const int length)
	{
		Impl identityMatrixImpl(length, length);
		for (size_t row = 0; row < identityMatrixImpl.height(); row++) {
			for (size_t col = 0; col < identityMatrixImpl.width(); col++) {
				identityMatrixImpl[row][col] = (row == col) ? 1.0 : 0.0;
			}
		}
		return identityMatrixImpl;
	}

	const Roww& Matrixx::Impl::operator[](const size_t row) const
	{
		int exceptNum = ExceptionHandlerr::checkRowIndex(row, mHeight);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		return mRows[row];
	}
	Roww& Matrixx::Impl::operator[](const size_t row)
	{
		return const_cast<Roww&>(static_cast<const Impl&>(*this)[row]);
	}

	const Roww& Matrixx::Impl::operator()(const int row) const
	{
		int exceptNum = ExceptionHandlerr::checkRowIndex(row, mHeight);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight, true);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
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
	Roww& Matrixx::Impl::operator()(const int row)
	{
		return const_cast<Roww&>(static_cast<const Impl&>(*this)(row));
	}

	const double& Matrixx::Impl::operator()(const int row, const int col) const
	{
		int exceptNum = ExceptionHandlerr::checkRowIndex(row, mHeight);
		exceptNum += ExceptionHandlerr::checkColumnIndex(col, mWidth);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight, true);
			ColumnIndexArgument colIndexArg(col, mWidth, true);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
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
	double& Matrixx::Impl::operator()(const int row, const int col)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this)(row, col));
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

	Matrixx::Impl Matrixx::Impl::operator+() const
	{
		return Impl(*this);
	}
	Matrixx::Impl Matrixx::Impl::operator-() const
	{
		Impl negativeMatrixImpl(static_cast<int>(mHeight), static_cast<int>(mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				negativeMatrixImpl[row][col] = convertNegativeZero(-mRows[row][col]);
			}
		}
		return negativeMatrixImpl;
	}

	void Matrixx::Impl::swap(Impl& rightMatrixImpl) noexcept
	{
		std::swap(mRows, rightMatrixImpl.mRows);
		std::swap(mHeight, rightMatrixImpl.mHeight);
		std::swap(mWidth, rightMatrixImpl.mWidth);
	}

	Matrixx::Impl& Matrixx::Impl::operator+=(const Impl& rightMatrixImpl)
	{
		int exceptNum = ExceptionHandlerr::checkHeight(mHeight, rightMatrixImpl.mHeight);
		exceptNum += ExceptionHandlerr::checkWidth(mWidth, rightMatrixImpl.mWidth);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightMatrixImpl.mHeight, rightMatrixImpl.mWidth);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultMatrixImpl(*this);
		for (size_t row = 0; row < mHeight; row++) {
			resultMatrixImpl[row] += rightMatrixImpl[row];
		}
		swap(resultMatrixImpl);
		return *this;
	}
	Matrixx::Impl& Matrixx::Impl::operator-=(const Impl& rightMatrixImpl)
	{
		int exceptNum = ExceptionHandlerr::checkHeight(mHeight, rightMatrixImpl.mHeight);
		exceptNum += ExceptionHandlerr::checkWidth(mWidth, rightMatrixImpl.mWidth);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightMatrixImpl.mHeight, rightMatrixImpl.mWidth);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultMatrixImpl(*this);
		for (size_t row = 0; row < mHeight; row++) {
			resultMatrixImpl[row] -= rightMatrixImpl[row];
		}
		swap(resultMatrixImpl);
		return *this;
	}
	Matrixx::Impl& Matrixx::Impl::operator*=(const double multiplier)
	{
		for (size_t row = 0; row < mHeight; row++) {
			mRows[row] *= multiplier;
		}
		return *this;
	}
	Matrixx::Impl& Matrixx::Impl::operator*=(const Impl& rightMatrixImpl)
	{
		int exceptNum = ExceptionHandlerr::checkJoinLength(mWidth, rightMatrixImpl.mHeight);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightMatrixImpl.mHeight, rightMatrixImpl.mWidth);
			OperationArgument operationArg('*', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultMatrixImpl(static_cast<int>(mHeight), static_cast<int>(rightMatrixImpl.mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < rightMatrixImpl.mWidth; col++) {
				double dotProduct = 0.0;
				for (size_t join = 0; join < mWidth; join++) {
					dotProduct += mRows[row][join] * rightMatrixImpl[join][col];
				}
				resultMatrixImpl[row][col] = convertNegativeZero(dotProduct);
			}
		}
		swap(resultMatrixImpl);
		return *this;
	}
	Matrixx::Impl& Matrixx::Impl::operator/=(const double divisor)
	{
		if (convertNegativeZero(divisor) == 0.0) {
			ExceptionHandlerr handler(ExceptionState::ArithmeticException,
				static_cast<int>(OperationState::DivideByZero));
			handler.handleException();
		}

		Impl resultMatrixImpl(*this);
		for (size_t row = 0; row < mHeight; row++) {
			resultMatrixImpl[row] /= divisor;
		}
		swap(resultMatrixImpl);
		return *this;
	}

	Matrixx::Impl& Matrixx::Impl::operator&=(const Impl& rightMatrixImpl)
	{
		int exceptNum = ExceptionHandlerr::checkHeight(mHeight, rightMatrixImpl.mHeight);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightMatrixImpl.mHeight, rightMatrixImpl.mWidth);
			OperationArgument operationArg('&', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl appendedMatrixImpl(static_cast<int>(mHeight), static_cast<int>(mWidth + rightMatrixImpl.mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrixImpl[row][col] = mRows[row][col];
			}
		}
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = mWidth; col < mWidth + rightMatrixImpl.mWidth; col++) {
				appendedMatrixImpl[row][col] = rightMatrixImpl[row][col - mWidth];
			}
		}
		swap(appendedMatrixImpl);
		return *this;
	}
	Matrixx::Impl& Matrixx::Impl::operator&=(const Vectorr::Impl& rightVectorImpl)
	{
		return (*this &= Impl(rightVectorImpl));
	}
	
	Matrixx::Impl& Matrixx::Impl::operator|=(const Impl& lowerMatrixImpl)
	{
		int exceptNum = ExceptionHandlerr::checkWidth(mWidth, lowerMatrixImpl.mWidth);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument upperLengthArg(mHeight, mWidth);
			LengthArgument lowerLengthArg(lowerMatrixImpl.mHeight, lowerMatrixImpl.mWidth);
			OperationArgument operationArg('|', upperLengthArg, lowerLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl appendedMatrixImpl(static_cast<int>(mHeight + lowerMatrixImpl.mHeight), static_cast<int>(mWidth));
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrixImpl[row][col] = mRows[row][col];
			}
		}
		for (size_t row = mHeight; row < mHeight + lowerMatrixImpl.mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				appendedMatrixImpl[row][col] = lowerMatrixImpl[row - mHeight][col];
			}
		}
		swap(appendedMatrixImpl);
		return *this;
	}
	Matrixx::Impl& Matrixx::Impl::operator|=(const Roww::Impl& lowerRowImpl)
	{
		return (*this |= Impl(lowerRowImpl));
	}

	Vectorr::Impl Matrixx::Impl::operator*(const Vectorr::Impl& rightVectorImpl) const
	{
		int exceptNum = ExceptionHandlerr::checkJoinLength(mWidth, rightVectorImpl.mHeight);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightVectorImpl.mHeight, 1);
			OperationArgument operationArg('*', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Vectorr::Impl resultVectorImpl(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			double dotProduct = 0.0;
			for (size_t join = 0; join < mWidth; join++) {
				dotProduct += mRows[row][join] * rightVectorImpl[join];
			}
			resultVectorImpl[row] = dotProduct;
		}
		return resultVectorImpl;
	}

	Matrixx::Impl Matrixx::Impl::operator+(const Impl& rightMatrixImpl) const
	{
		Impl resultMatrixImpl(*this);
		resultMatrixImpl += rightMatrixImpl;
		return resultMatrixImpl;
	}
	Matrixx::Impl Matrixx::Impl::operator-(const Impl& rightMatrixImpl) const
	{
		Impl resultMatrixImpl(*this);
		resultMatrixImpl -= rightMatrixImpl;
		return resultMatrixImpl;
	}
	Matrixx::Impl Matrixx::Impl::operator*(const double multiplier) const
	{
		Impl resultMatrixImpl(*this);
		resultMatrixImpl *= multiplier;
		return resultMatrixImpl;
	}
	Matrixx::Impl Matrixx::Impl::operator*(const Impl& rightMatrixImpl) const
	{
		Impl resultMatrixImpl(*this);
		resultMatrixImpl *= rightMatrixImpl;
		return resultMatrixImpl;
	}
	Matrixx::Impl Matrixx::Impl::operator/(const double divisor) const
	{
		Impl resultMatrixImpl(*this);
		resultMatrixImpl /= divisor;
		return resultMatrixImpl;
	}

	Matrixx::Impl Matrixx::Impl::operator&(const Impl& rightMatrixImpl) const
	{
		Impl resultMatrixImpl(*this);
		resultMatrixImpl &= rightMatrixImpl;
		return resultMatrixImpl;
	}
	Matrixx::Impl Matrixx::Impl::operator&(const Vectorr::Impl& rightVectorImpl) const
	{
		Impl resultMatrixImpl(*this);
		resultMatrixImpl &= rightVectorImpl;
		return resultMatrixImpl;
	}

	Matrixx::Impl Matrixx::Impl::operator|(const Impl& lowerMatrixImpl) const
	{
		Impl resultMatrixImpl(*this);
		resultMatrixImpl |= lowerMatrixImpl;
		return resultMatrixImpl;
	}
	Matrixx::Impl Matrixx::Impl::operator|(const Roww::Impl& lowerRowImpl) const
	{
		Impl resultMatrixImpl(*this);
		resultMatrixImpl |= lowerRowImpl;
		return resultMatrixImpl;
	}

	bool Matrixx::Impl::operator==(const Impl& rightMatrixImpl) const
	{
		if (mHeight != rightMatrixImpl.mHeight ||
			mWidth != rightMatrixImpl.mWidth) {
			return false;
		}
		for (size_t row = 0; row < mHeight; row++) {
			if (mRows[row] != rightMatrixImpl[row]) {
				return false;
			}
		}
		return true;
	}
	bool Matrixx::Impl::operator!=(const Impl& rightMatrixImpl) const
	{
		return !(*this == rightMatrixImpl);
	}

	const Roww::Impl Matrixx::Impl::getRow(const int row) const
	{
		int exceptNum = ExceptionHandlerr::checkRowIndex(row, mHeight);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight, true);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		if (row >= 0) {
			return *(mRows[row].impl);
		}
		else {
			return *(mRows[static_cast<size_t>(static_cast<int>(mHeight) + row)].impl);
		}
	}
	const Vectorr::Impl Matrixx::Impl::getColumn(const int col) const
	{
		int exceptNum = ExceptionHandlerr::checkColumnIndex(col, mWidth);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(col, mWidth, true);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		Vectorr::Impl copyVectorImpl(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			copyVectorImpl[row] = mRows[row](col);
		}
		return copyVectorImpl;
	}

	const size_t Matrixx::Impl::height() const
	{
		return mHeight;
	}
	const size_t Matrixx::Impl::width() const
	{
		return mWidth;
	}
	const size_t Matrixx::Impl::size() const
	{
		return mHeight * mWidth;
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

	const double& Roww::Impl::operator[](const size_t col) const
	{
		int exceptNum = ExceptionHandlerr::checkColumnIndex(col, mWidth);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			ColumnIndexArgument colIndexArg(col, mWidth);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		return mEntries[col];
	}
	double& Roww::Impl::operator[](const size_t col)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this)[col]);
	}
	const double& Roww::Impl::operator()(const int col) const
	{
		int exceptNum = ExceptionHandlerr::checkColumnIndex(col, mWidth);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			ColumnIndexArgument colIndexArg(col, mWidth, true);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		if (col >= 0) {
			return mEntries[col];
		}
		else {
			return mEntries[static_cast<size_t>(static_cast<int>(mWidth) + col)];
		}
	}
	double& Roww::Impl::operator()(const int col)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this)(col));
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

	Roww::Impl Roww::Impl::operator+() const
	{
		return Impl(*this);
	}
	Roww::Impl Roww::Impl::operator-() const
	{
		Impl negativeRowImpl(static_cast<int>(mWidth));
		for (size_t col = 0; col < mWidth; col++) {
			negativeRowImpl[col] = convertNegativeZero(-mEntries[col]);
		}
		return negativeRowImpl;
	}

	void Roww::Impl::swap(Impl& rightRowImpl) noexcept
	{
		std::swap(mEntries, rightRowImpl.mEntries);
		std::swap(mWidth, rightRowImpl.mWidth);
	}

	Roww::Impl& Roww::Impl::operator+=(const Impl& rightRowImpl)
	{
		int exceptNum = ExceptionHandlerr::checkWidth(mWidth, rightRowImpl.mWidth);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(1, mWidth);
			LengthArgument rightLengthArg(1, rightRowImpl.mWidth);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultRowImpl(static_cast<int>(mWidth));
		for (size_t col = 0; col < mWidth; col++) {
			resultRowImpl[col] = convertNegativeZero(mEntries[col] + rightRowImpl[col]);
		}
		swap(resultRowImpl);
		return *this;
	}
	Roww::Impl& Roww::Impl::operator-=(const Impl& rightRowImpl)
	{
		int exceptNum = ExceptionHandlerr::checkWidth(mWidth, rightRowImpl.mWidth);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(1, mWidth);
			LengthArgument rightLengthArg(1, rightRowImpl.mWidth);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultRowImpl(static_cast<int>(mWidth));
		for (size_t col = 0; col < mWidth; col++) {
			resultRowImpl[col] = convertNegativeZero(mEntries[col] - rightRowImpl[col]);
		}
		swap(resultRowImpl);
		return *this;
	}
	Roww::Impl& Roww::Impl::operator*=(const double multiplier)
	{
		for (size_t col = 0; col < mWidth; col++) {
			mEntries[col] = convertNegativeZero(multiplier * mEntries[col]);
		}
		return *this;
	}
	Roww::Impl& Roww::Impl::operator/=(const double divisor)
	{
		if (convertNegativeZero(divisor) == 0.0) {
			ExceptionHandlerr handler(ExceptionState::ArithmeticException,
				static_cast<int>(OperationState::DivideByZero));
			handler.handleException();
		}

		Impl resultRowImpl(static_cast<int>(mWidth));
		for (size_t col = 0; col < mWidth; col++) {
			resultRowImpl[col] = convertNegativeZero(mEntries[col] / divisor);
		}
		swap(resultRowImpl);
		return *this;
	}
	Roww::Impl& Roww::Impl::operator&=(const Impl& rightRowImpl)
	{
		Impl appendedRowImpl(static_cast<int>(mWidth + rightRowImpl.mWidth));
		for (size_t col = 0; col < mWidth; col++) {
			appendedRowImpl[col] = mEntries[col];
		}
		for (size_t col = mWidth; col < mWidth + rightRowImpl.mWidth; col++) {
			appendedRowImpl[col] = rightRowImpl[col - mWidth];
		}
		swap(appendedRowImpl);
		return *this;
	}

	Roww::Impl Roww::Impl::operator+(const Impl& rightRowImpl) const
	{
		Impl resultRowImpl(*this);
		resultRowImpl += rightRowImpl;
		return resultRowImpl;
	}
	Roww::Impl Roww::Impl::operator-(const Impl& rightRowImpl) const
	{
		Impl resultRowImpl(*this);
		resultRowImpl -= rightRowImpl;
		return resultRowImpl;
	}
	Roww::Impl Roww::Impl::operator*(const double multiplier) const
	{
		Impl resultRowImpl(*this);
		resultRowImpl *= multiplier;
		return resultRowImpl;
	}
	Roww::Impl Roww::Impl::operator/(const double divisor) const
	{
		Impl resultRowImpl(*this);
		resultRowImpl /= divisor;
		return resultRowImpl;
	}

	Roww::Impl Roww::Impl::operator&(const Impl& rightRowImpl) const
	{
		Impl resultRowImpl(*this);
		resultRowImpl &= rightRowImpl;
		return resultRowImpl;
	}

	Matrixx::Impl Roww::Impl::operator|(const Matrixx::Impl& lowerMatrixImpl) const
	{
		Matrixx::Impl resultMatrixImpl(*this);
		resultMatrixImpl |= lowerMatrixImpl;
		return resultMatrixImpl;
	}
	Matrixx::Impl Roww::Impl::operator|(const Roww::Impl& lowerRowImpl) const
	{
		Matrixx::Impl resultMatrixImpl(*this);
		resultMatrixImpl |= lowerRowImpl;
		return resultMatrixImpl;
	}

	bool Roww::Impl::operator==(const Impl& rightRowImpl) const
	{
		if (mWidth != rightRowImpl.mWidth) {
			return false;
		}
		for (size_t col = 0; col < mWidth; col++) {
			if (mEntries[col] != rightRowImpl[col]) {
				return false;
			}
		}
		return true;
	}
	bool Roww::Impl::operator!=(const Impl& rightRowImpl) const
	{
		return !(*this == rightRowImpl);
	}

	const size_t Roww::Impl::width() const
	{
		return mWidth;
	}
	const size_t Roww::Impl::size() const
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

	const double& Vectorr::Impl::operator[](const size_t row) const
	{
		int exceptNum = ExceptionHandlerr::checkRowIndex(row, mHeight);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		return mEntries[row];
	}
	double& Vectorr::Impl::operator[](const size_t row)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this)[row]);
	}
	const double& Vectorr::Impl::operator()(const int row) const
	{
		int exceptNum = ExceptionHandlerr::checkRowIndex(row, mHeight);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mHeight, true);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		if (row >= 0) {
			return mEntries[row];
		}
		else {
			return mEntries[static_cast<size_t>(static_cast<int>(mHeight) + row)];
		}
	}
	double& Vectorr::Impl::operator()(const int row)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this)(row));
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

	Vectorr::Impl Vectorr::Impl::operator+() const
	{
		return Impl(*this);
	}
	Vectorr::Impl Vectorr::Impl::operator-() const
	{
		Impl negativeVectorImpl(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			negativeVectorImpl[row] = convertNegativeZero(-mEntries[row]);
		}
		return negativeVectorImpl;
	}

	void Vectorr::Impl::swap(Impl& rightVectorImpl) noexcept
	{
		std::swap(mEntries, rightVectorImpl.mEntries);
		std::swap(mHeight, rightVectorImpl.mHeight);
	}

	Vectorr::Impl& Vectorr::Impl::operator+=(const Vectorr::Impl& rightVectorImpl)
	{
		int exceptNum = ExceptionHandlerr::checkHeight(mHeight, rightVectorImpl.mHeight);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, 1);
			LengthArgument rightLengthArg(rightVectorImpl.mHeight, 1);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultVectorImpl(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			resultVectorImpl[row] = convertNegativeZero(mEntries[row] + rightVectorImpl[row]);
		}
		swap(resultVectorImpl);
		return *this;
	}
	Vectorr::Impl& Vectorr::Impl::operator-=(const Vectorr::Impl& rightVectorImpl)\
	{
		int exceptNum = ExceptionHandlerr::checkHeight(mHeight, rightVectorImpl.mHeight);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, 1);
			LengthArgument rightLengthArg(rightVectorImpl.mHeight, 1);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultVectorImpl(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			resultVectorImpl[row] = convertNegativeZero(mEntries[row] - rightVectorImpl[row]);
		}
		swap(resultVectorImpl);
		return *this;
	}
	Vectorr::Impl& Vectorr::Impl::operator*=(const double multiplier)
	{
		for (size_t row = 0; row < mHeight; row++) {
			mEntries[row] = convertNegativeZero(multiplier * mEntries[row]);
		}
		return *this;
	}
	Vectorr::Impl& Vectorr::Impl::operator/=(const double divisor)
	{
		if (convertNegativeZero(divisor) == 0.0) {
			ExceptionHandlerr handler(ExceptionState::ArithmeticException,
				static_cast<int>(OperationState::DivideByZero));
			handler.handleException();
		}

		Impl resultVectorImpl(static_cast<int>(mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			resultVectorImpl[row] = convertNegativeZero(mEntries[row] / divisor);
		}
		swap(resultVectorImpl);
		return *this;
	}

	Vectorr::Impl& Vectorr::Impl::operator|=(const Vectorr::Impl& lowerVectorImpl)
	{
		Impl appendedVectorImpl(static_cast<int>(mHeight + lowerVectorImpl.mHeight));
		for (size_t row = 0; row < mHeight; row++) {
			appendedVectorImpl[row] = mEntries[row];
		}
		for (size_t row = mHeight; row < mHeight + lowerVectorImpl.mHeight; row++) {
			appendedVectorImpl[row] = lowerVectorImpl[row - mHeight];
		}
		swap(appendedVectorImpl);
		return *this;
	}

	Vectorr::Impl Vectorr::Impl::operator+(const Vectorr::Impl& rightVectorImpl) const
	{
		Impl resultVectorImpl(*this);
		resultVectorImpl += rightVectorImpl;
		return resultVectorImpl;
	}
	Vectorr::Impl Vectorr::Impl::operator-(const Vectorr::Impl& rightVectorImpl) const
	{
		Impl resultVectorImpl(*this);
		resultVectorImpl -= rightVectorImpl;
		return resultVectorImpl;
	}
	Vectorr::Impl Vectorr::Impl::operator*(const double multiplier) const
	{
		Impl resultVectorImpl(*this);
		resultVectorImpl *= multiplier;
		return resultVectorImpl;
	}
	Vectorr::Impl Vectorr::Impl::operator/(const double divisor) const
	{
		Impl resultVectorImpl(*this);
		resultVectorImpl /= divisor;
		return resultVectorImpl;
	}

	Matrixx::Impl Vectorr::Impl::operator&(const Matrixx::Impl& rightMatrixImpl) const
	{
		Matrixx::Impl resultMatrixImpl(*this);
		resultMatrixImpl &= rightMatrixImpl;
		return resultMatrixImpl;
	}
	Matrixx::Impl Vectorr::Impl::operator&(const Impl& rightVectorImpl) const
	{
		Matrixx::Impl resultMatrixImpl(*this);
		resultMatrixImpl &= rightVectorImpl;
		return resultMatrixImpl;
	}

	Vectorr::Impl Vectorr::Impl::operator|(const Vectorr::Impl& lowerVectorImpl) const
	{
		Impl resultVectorImpl(*this);
		resultVectorImpl |= lowerVectorImpl;
		return resultVectorImpl;
	}

	bool Vectorr::Impl::operator==(const Vectorr::Impl& rightVectorImpl) const
	{
		if (mHeight != rightVectorImpl.mHeight) {
			return false;
		}
		for (size_t row = 0; row < mHeight; row++) {
			if (mEntries[row] != rightVectorImpl[row]) {
				return false;
			}
		}
		return true;
	}
	bool Vectorr::Impl::operator!=(const Vectorr::Impl& rightVectorImpl) const
	{
		return !(*this == rightVectorImpl);
	}

	const size_t Vectorr::Impl::height() const
	{
		return mHeight;
	}
	const size_t Vectorr::Impl::size() const
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