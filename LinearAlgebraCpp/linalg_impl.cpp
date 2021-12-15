#include "linalg_impl.h"

namespace linalg {

	Tensorr::Impl::Impl(const int size)
		: mSize(size)
	{
	}

	const size_t Tensorr::Impl::size() const
	{
		return mSize;
	}
	void Tensorr::Impl::size(size_t size)
	{
		mSize = size;
	}

	const double Tensorr::Impl::epsilonTest(const double value)
	{
		constexpr double epsilon = std::numeric_limits<double>::epsilon();
		return (value < epsilon && value > -epsilon) ? 0.0 : value;
	}





	Matrixx::Impl::Impl(const size_t height, const size_t width)
	{
		init(height, width);
	}
	Matrixx::Impl::Impl(const Roww::Impl& copyRowImpl)
		: Impl(1, copyRowImpl.mSize)
	{
		for (size_t col = 0; col < mHeight; col++) {
			mRows[0][col] = copyRowImpl[col];
		}
	}
	Matrixx::Impl::Impl(const Vectorr::Impl& copyVectorImpl)
		: Impl(copyVectorImpl.mSize, 1)
	{
		for (size_t row = 0; row < mHeight; row++) {
			mRows[row][0] = copyVectorImpl[row];
		}
	}
	void Matrixx::Impl::init(const size_t height, const size_t width)
	{
		int exceptNum = ExceptionHandlerr::checkValidHeight(height);
		exceptNum += ExceptionHandlerr::checkValidWidth(width);
		if (exceptNum > static_cast<int>(LengthState::NoExcept)) {
			LengthArgument lengthArg(height, width);
			ExceptionHandlerr handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}

		Tensorr::Impl::size(static_cast<size_t>(height * width));

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

		for (int row = static_cast<int>(mHeight - 1); row >= 0; row--) {
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

		Impl blockMatrixImpl(blockHeight, blockWidth);
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

		const size_t length = mHeight;
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
		Impl transposedMatrixImpl(mWidth, mHeight);
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				transposedMatrixImpl[col][row] = mRows[row][col];
			}
		}
		return transposedMatrixImpl;
	}

	Matrixx::Impl Matrixx::Impl::identity(const size_t length)
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
			mRows[sequence / mWidth][sequence % mWidth] = epsilonTest(value);
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
		Impl negativeMatrixImpl(mHeight, mWidth);
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < mWidth; col++) {
				negativeMatrixImpl[row][col] = epsilonTest(-mRows[row][col]);
			}
		}
		return negativeMatrixImpl;
	}

	void Matrixx::Impl::swap(Impl& rightMatrixImpl) noexcept
	{
		std::swap(mRows, rightMatrixImpl.mRows);

		std::swap(mSize, rightMatrixImpl.mSize);

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

		Impl resultMatrixImpl(mHeight, rightMatrixImpl.mWidth);
		for (size_t row = 0; row < mHeight; row++) {
			for (size_t col = 0; col < rightMatrixImpl.mWidth; col++) {
				double dotProduct = 0.0;
				for (size_t join = 0; join < mWidth; join++) {
					dotProduct += mRows[row][join] * rightMatrixImpl[join][col];
				}
				resultMatrixImpl[row][col] = epsilonTest(dotProduct);
			}
		}
		swap(resultMatrixImpl);
		return *this;
	}
	Matrixx::Impl& Matrixx::Impl::operator/=(const double divisor)
	{
		if (epsilonTest(divisor) == 0.0) {
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

		Impl appendedMatrixImpl(mHeight, mWidth + rightMatrixImpl.mWidth);
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
		return *this &= Impl(rightVectorImpl);
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

		Impl appendedMatrixImpl(mHeight + lowerMatrixImpl.mHeight, mWidth);
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
		return *this |= Impl(lowerRowImpl);
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

	Vectorr::Impl Matrixx::Impl::operator*(const Vectorr::Impl& rightVectorImpl) const
	{
		int exceptNum = ExceptionHandlerr::checkJoinLength(mWidth, rightVectorImpl.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mHeight, mWidth);
			LengthArgument rightLengthArg(rightVectorImpl.mSize, 1);
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

		Vectorr::Impl copyVectorImpl(mHeight);
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

	const std::string Matrixx::Impl::str() const
	{
		std::string matrixString = "(" + std::to_string(mHeight) + " x " + std::to_string(mWidth) + " Matrix)\n";
		for (size_t row = 0; row < mHeight; row++) {
			matrixString += mRows[row].str();
		}
		return matrixString;
	}





	Roww::Impl::Impl(const size_t size)
	{
		init(size);
	}
	void Roww::Impl::init(const size_t size)
	{
		int exceptNum = ExceptionHandlerr::checkValidWidth(size);
		if (exceptNum > static_cast<int>(LengthState::NoExcept)) {
			LengthArgument lengthArg(1, size);
			ExceptionHandlerr handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}
		
		Tensorr::Impl::size(size);

		mEntries.clear();
		mEntries.resize(mSize, 0.0);
	}

	const double& Roww::Impl::operator[](const size_t col) const
	{
		int exceptNum = ExceptionHandlerr::checkColumnIndex(col, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			ColumnIndexArgument colIndexArg(col, mSize);
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
		int exceptNum = ExceptionHandlerr::checkColumnIndex(col, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			ColumnIndexArgument colIndexArg(col, mSize, true);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
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
	double& Roww::Impl::operator()(const int col)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this)(col));
	}

	void Roww::Impl::allocate(const size_t sequence, const double value)
	{
		if (sequence < mSize) {
			mEntries[sequence] = epsilonTest(value);
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
		Impl negativeRowImpl(mSize);
		for (size_t col = 0; col < mSize; col++) {
			negativeRowImpl[col] = epsilonTest(-mEntries[col]);
		}
		return negativeRowImpl;
	}

	void Roww::Impl::swap(Impl& rightRowImpl) noexcept
	{
		std::swap(mEntries, rightRowImpl.mEntries);

		std::swap(mSize, rightRowImpl.mSize);
	}

	Roww::Impl& Roww::Impl::operator+=(const Impl& rightRowImpl)
	{
		int exceptNum = ExceptionHandlerr::checkWidth(mSize, rightRowImpl.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(1, mSize);
			LengthArgument rightLengthArg(1, rightRowImpl.mSize);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultRowImpl(mSize);
		for (size_t col = 0; col < mSize; col++) {
			resultRowImpl[col] = epsilonTest(mEntries[col] + rightRowImpl[col]);
		}
		swap(resultRowImpl);
		return *this;
	}
	Roww::Impl& Roww::Impl::operator-=(const Impl& rightRowImpl)
	{
		int exceptNum = ExceptionHandlerr::checkWidth(mSize, rightRowImpl.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(1, mSize);
			LengthArgument rightLengthArg(1, rightRowImpl.mSize);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultRowImpl(mSize);
		for (size_t col = 0; col < mSize; col++) {
			resultRowImpl[col] = epsilonTest(mEntries[col] - rightRowImpl[col]);
		}
		swap(resultRowImpl);
		return *this;
	}
	Roww::Impl& Roww::Impl::operator*=(const double multiplier)
	{
		for (size_t col = 0; col < mSize; col++) {
			mEntries[col] = epsilonTest(multiplier * mEntries[col]);
		}
		return *this;
	}
	Roww::Impl& Roww::Impl::operator/=(const double divisor)
	{
		if (epsilonTest(divisor) == 0.0) {
			ExceptionHandlerr handler(ExceptionState::ArithmeticException,
				static_cast<int>(OperationState::DivideByZero));
			handler.handleException();
		}

		Impl resultRowImpl(mSize);
		for (size_t col = 0; col < mSize; col++) {
			resultRowImpl[col] = epsilonTest(mEntries[col] / divisor);
		}
		swap(resultRowImpl);
		return *this;
	}

	Roww::Impl& Roww::Impl::operator&=(const Impl& rightRowImpl)
	{
		Impl appendedRowImpl(mSize + rightRowImpl.mSize);
		for (size_t col = 0; col < mSize; col++) {
			appendedRowImpl[col] = mEntries[col];
		}
		for (size_t col = mSize; col < mSize + rightRowImpl.mSize; col++) {
			appendedRowImpl[col] = rightRowImpl[col - mSize];
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
		if (mSize != rightRowImpl.mSize) {
			return false;
		}
		for (size_t col = 0; col < mSize; col++) {
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

	const std::string Roww::Impl::str() const
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





	Vectorr::Impl::Impl(const size_t size)
	{
		init(size);
	}
	void Vectorr::Impl::init(const size_t size)
	{
		int exceptNum = ExceptionHandlerr::checkValidHeight(size);
		if (exceptNum > static_cast<int>(LengthState::NoExcept)) {
			LengthArgument lengthArg(size, 1);
			ExceptionHandlerr handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}
		
		Tensorr::Impl::size(size);

		mEntries.clear();
		mEntries.resize(mSize, 0.0);
	}

	const double& Vectorr::Impl::operator[](const size_t row) const
	{
		int exceptNum = ExceptionHandlerr::checkRowIndex(row, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mSize);
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
		int exceptNum = ExceptionHandlerr::checkRowIndex(row, mSize);
		if (exceptNum > static_cast<int>(IndexState::NoExcept)) {
			RowIndexArgument rowIndexArg(row, mSize, true);
			ExceptionHandlerr handler(ExceptionState::OutOfRange, exceptNum);
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
	double& Vectorr::Impl::operator()(const int row)
	{
		return const_cast<double&>(static_cast<const Impl&>(*this)(row));
	}

	void Vectorr::Impl::allocate(const size_t sequence, const double value)
	{
		if (sequence < mSize) {
			mEntries[sequence] = epsilonTest(value);
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
		Impl negativeVectorImpl(mSize);
		for (size_t row = 0; row < mSize; row++) {
			negativeVectorImpl[row] = epsilonTest(-mEntries[row]);
		}
		return negativeVectorImpl;
	}

	void Vectorr::Impl::swap(Impl& rightVectorImpl) noexcept
	{
		std::swap(mEntries, rightVectorImpl.mEntries);

		std::swap(mSize, rightVectorImpl.mSize);
	}

	Vectorr::Impl& Vectorr::Impl::operator+=(const Vectorr::Impl& rightVectorImpl)
	{
		int exceptNum = ExceptionHandlerr::checkHeight(mSize, rightVectorImpl.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mSize, 1);
			LengthArgument rightLengthArg(rightVectorImpl.mSize, 1);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultVectorImpl(mSize);
		for (size_t row = 0; row < mSize; row++) {
			resultVectorImpl[row] = epsilonTest(mEntries[row] + rightVectorImpl[row]);
		}
		swap(resultVectorImpl);
		return *this;
	}
	Vectorr::Impl& Vectorr::Impl::operator-=(const Vectorr::Impl& rightVectorImpl)\
	{
		int exceptNum = ExceptionHandlerr::checkHeight(mSize, rightVectorImpl.mSize);
		if (exceptNum > static_cast<int>(OperationState::NoExcept)) {
			LengthArgument leftLengthArg(mSize, 1);
			LengthArgument rightLengthArg(rightVectorImpl.mSize, 1);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandlerr handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}

		Impl resultVectorImpl(mSize);
		for (size_t row = 0; row < mSize; row++) {
			resultVectorImpl[row] = epsilonTest(mEntries[row] - rightVectorImpl[row]);
		}
		swap(resultVectorImpl);
		return *this;
	}
	Vectorr::Impl& Vectorr::Impl::operator*=(const double multiplier)
	{
		for (size_t row = 0; row < mSize; row++) {
			mEntries[row] = epsilonTest(multiplier * mEntries[row]);
		}
		return *this;
	}
	Vectorr::Impl& Vectorr::Impl::operator/=(const double divisor)
	{
		if (epsilonTest(divisor) == 0.0) {
			ExceptionHandlerr handler(ExceptionState::ArithmeticException,
				static_cast<int>(OperationState::DivideByZero));
			handler.handleException();
		}

		Impl resultVectorImpl(mSize);
		for (size_t row = 0; row < mSize; row++) {
			resultVectorImpl[row] = epsilonTest(mEntries[row] / divisor);
		}
		swap(resultVectorImpl);
		return *this;
	}

	Vectorr::Impl& Vectorr::Impl::operator|=(const Vectorr::Impl& lowerVectorImpl)
	{
		Impl appendedVectorImpl(mSize + lowerVectorImpl.mSize);
		for (size_t row = 0; row < mSize; row++) {
			appendedVectorImpl[row] = mEntries[row];
		}
		for (size_t row = mSize; row < mSize + lowerVectorImpl.mSize; row++) {
			appendedVectorImpl[row] = lowerVectorImpl[row - mSize];
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
		if (mSize != rightVectorImpl.mSize) {
			return false;
		}
		for (size_t row = 0; row < mSize; row++) {
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

	const std::string Vectorr::Impl::str() const
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