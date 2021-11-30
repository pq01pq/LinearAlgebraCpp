#include "linalg.h"

namespace linalg {
	Allocator::Allocator(Allocatable& target, const int sequence)
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





	Matrixx::Matrixx(const int height, const int width)
	{
		check::LengthInfo lengthInfo{ height, width };
		int errorNumber = check::checkHeight(height);
		errorNumber += check::checkWidth(width);
		check::handleLengthError(errorNumber, lengthInfo);

		this->height = height;
		this->width = width;
		rows = new Roww[height];
		for (int row = 0; row < height; row++) {
			rows[row].init(width);
		}
	}
	Matrixx::Matrixx(const Matrixx& copyMatrix)
		: Matrixx(copyMatrix.height, copyMatrix.width)
	{
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				rows[row][col] = copyMatrix[row][col];
			}
		}
	}
	Matrixx::Matrixx(const Roww& copyRow)
		: Matrixx(1, copyRow.width)
	{
		for (int col = 0; col < height; col++) {
			rows[0][col] = copyRow[col];
		}
	}
	Matrixx::Matrixx(const Vectorr& copyVector)
		: Matrixx(copyVector.height, 1)
	{
		for (int row = 0; row < height; row++) {
			rows[row][0] = copyVector[row];
		}
	}
	Matrixx::~Matrixx()
	{
		delete[] rows;
	}

	Matrixx Matrixx::block(const int beginRow, const int beginCol, const int blockHeight, const int blockWidth) const
	{
		check::LengthInfo lengthInfo{ height, width };
		check::IndexInfo beginIndexInfo{ beginRow, beginCol };
		check::IndexInfo endIndexInfo{ beginRow + blockHeight - 1, beginCol + blockWidth - 1 };

		int errorNumber = check::checkRowIndex(beginRow, height);
		errorNumber += check::checkColumnIndex(beginCol, width);
		check::handleOutOfRange(errorNumber, beginIndexInfo);

		errorNumber = check::checkRowIndex(beginRow + blockHeight - 1, height);
		errorNumber += check::checkColumnIndex(beginCol + blockWidth - 1, width);
		check::handleOutOfRange(errorNumber, endIndexInfo);

		Matrixx blockMatrix(blockHeight, blockWidth);
		for (int row = 0; row < blockHeight; row++) {
			for (int col = 0; col < blockWidth; col++) {
				blockMatrix[row][col] = rows[beginRow + row][beginCol + col];
			}
		}
		return blockMatrix;
	}

	Roww& Matrixx::operator[](const int row)
	{
		check::LengthInfo lengthInfo{ height, width };
		check::IndexInfo indexInfo{ row, 0, lengthInfo };
		int errorNumber = check::checkRowIndex(row, height);
		check::handleOutOfRange(errorNumber, indexInfo);

		return rows[row];
	}
	const Roww& Matrixx::operator[](const int row) const
	{
		check::LengthInfo lengthInfo{ height, width };
		check::IndexInfo indexInfo{ row, 0, lengthInfo };
		int errorNumber = check::checkRowIndex(row, height);
		check::handleOutOfRange(errorNumber, indexInfo);

		return rows[row];
	}

	double& Matrixx::operator()(const int row, const int col) const
	{
		check::LengthInfo lengthInfo{ height, width };
		check::IndexInfo indexInfo{ row, col, lengthInfo };
		int errorNumber = check::checkRowIndex(row, height);
		check::handleOutOfRange(errorNumber, indexInfo);

		return rows[row][col];
	}

	Allocator& Matrixx::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Matrixx::allocate(const int sequence, const double value)
	{
		if (sequence < height * width) {
			rows[sequence / width][sequence % width] = check::preventNegativeZero(value);
		}
	}

	Matrixx Matrixx::operator+() const
	{
		return Matrixx(*this);
	}
	Matrixx Matrixx::operator-() const
	{
		Matrixx negativeMatrix(*this);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				negativeMatrix[row][col] = check::preventNegativeZero(-rows[row][col]);
			}
		}
		return negativeMatrix;
	}

	Matrixx& Matrixx::operator=(const Matrixx& rightMatrix)
	{
		if (this == &rightMatrix) {
			return *this;
		}

		Matrixx copyMatrix(rightMatrix);
		swap(*this, copyMatrix);
		return *this;
	}
	void swap(Matrixx& leftMatrix, Matrixx& rightMatrix) noexcept
	{
		std::swap(leftMatrix.rows, rightMatrix.rows);
		std::swap(leftMatrix.height, rightMatrix.height);
		std::swap(leftMatrix.width, rightMatrix.width);
	}

	Matrixx& Matrixx::operator+=(const Matrixx& rightMatrix)
	{
		check::LengthInfo leftLengthInfo{ height, width };
		check::LengthInfo rightLengthInfo{ rightMatrix.height, rightMatrix.width };
		check::OperationInfo operationInfo{ '+', leftLengthInfo, rightLengthInfo };
		int errorNumber = check::checkHeight(height, rightMatrix.height);
		errorNumber += check::checkWidth(width, rightMatrix.width);
		check::handleLogicError(errorNumber, operationInfo);
		
		for (int row = 0; row < height; row++) {
			rows[row] += rightMatrix[row];
		}
		return *this;
	}
	Matrixx& Matrixx::operator-=(const Matrixx& rightMatrix)
	{
		check::LengthInfo leftLengthInfo{ height, width };
		check::LengthInfo rightLengthInfo{ rightMatrix.height, rightMatrix.width };
		check::OperationInfo operationInfo{ '-', leftLengthInfo, rightLengthInfo };
		int errorNumber = check::checkHeight(height, rightMatrix.height);
		errorNumber += check::checkWidth(width, rightMatrix.width);
		check::handleLogicError(errorNumber, operationInfo);
		
		for (int row = 0; row < height; row++) {
			rows[row] -= rightMatrix[row];
		}
		return *this;
	}
	Matrixx& Matrixx::operator*=(const double multiplier)
	{
		for (int row = 0; row < height; row++) {
			rows[row] *= multiplier;
		}
		return *this;
	}

	Matrixx& Matrixx::operator*=(const Matrixx& rightMatrix)
	{
		check::LengthInfo leftLengthInfo{ height, width };
		check::LengthInfo rightLengthInfo{ rightMatrix.height, rightMatrix.width };
		check::OperationInfo operationInfo{ '*', leftLengthInfo, rightLengthInfo };
		int errorNumber = check::checkJoinLength(width, rightMatrix.height);
		check::handleLogicError(errorNumber, operationInfo);
		
		Matrixx resultMatrix(height, rightMatrix.width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < rightMatrix.width; col++) {
				double product = 0.0;
				for (int join = 0; join < width; join++) {
					product += rows[row][join] * rightMatrix[join][col];
				}
				resultMatrix[row][col] = check::preventNegativeZero(product);
			}
		}
		linalg::swap(*this, resultMatrix);
		return *this;
	}

	Matrixx& Matrixx::operator&=(const Matrixx& rightMatrix)
	{
		check::LengthInfo leftLengthInfo{ height, width };
		check::LengthInfo rightLengthInfo{ rightMatrix.height, rightMatrix.width };
		check::OperationInfo operationInfo{ '&', leftLengthInfo, rightLengthInfo };
		int errorNumber = check::checkHeight(height, rightMatrix.height);
		check::handleLogicError(errorNumber, operationInfo);
		
		Matrixx appendedMatrix(height, width + rightMatrix.width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				appendedMatrix[row][col] = rows[row][col];
			}
		}
		for (int row = 0; row < height; row++) {
			for (int col = width; col < width + rightMatrix.width; col++) {
				appendedMatrix[row][col] = rightMatrix[row][col - width];
			}
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator&=(const Vectorr& rightVector)
	{
		check::LengthInfo leftLengthInfo{ height, width };
		check::LengthInfo rightLengthInfo{ rightVector.height, 1 };
		check::OperationInfo operationInfo{ '&', leftLengthInfo, rightLengthInfo };
		int errorNumber = check::checkHeight(height, rightVector.height);
		check::handleLogicError(errorNumber, operationInfo);
		
		Matrixx appendedMatrix(height, width + 1);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				appendedMatrix[row][col] = rows[row][col];
			}
		}
		for (int row = 0; row < height; row++) {
			appendedMatrix[row][width] = rightVector[row];
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator|=(const Matrixx& lowerMatrix)
	{
		check::LengthInfo upperLengthInfo{ height, width };
		check::LengthInfo lowerLengthInfo{ lowerMatrix.height, lowerMatrix.width };
		check::OperationInfo operationInfo{ '|', upperLengthInfo, lowerLengthInfo };
		int errorNumber = check::checkWidth(width, lowerMatrix.width);
		check::handleLogicError(errorNumber, operationInfo);
		
		Matrixx appendedMatrix(height + lowerMatrix.height, width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				appendedMatrix[row][col] = rows[row][col];
			}
		}
		for (int row = height; row < height + lowerMatrix.height; row++) {
			for (int col = 0; col < width; col++) {
				appendedMatrix[row][col] = lowerMatrix[row - height][col];
			}
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator|=(const Roww& lowerRow)
	{
		check::LengthInfo upperLengthInfo{ height, width };
		check::LengthInfo lowerLengthInfo{ 1, lowerRow.width };
		check::OperationInfo operationInfo{ '|', upperLengthInfo, lowerLengthInfo };
		int errorNumber = check::checkWidth(width, lowerRow.width);
		check::handleLogicError(errorNumber, operationInfo);
		
		Matrixx appendedMatrix(height + 1, width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				appendedMatrix[row][col] = rows[row][col];
			}
		}
		for (int col = 0; col < width; col++) {
			appendedMatrix[height][col] = lowerRow[col];
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
	std::ostream& operator<<(std::ostream& outputStream, const Matrixx& outputMatrix)
	{
		outputStream << outputMatrix.str();
		return outputStream;
	}
	

	const int Matrixx::getHeight() const
	{
		return height;
	}
	const int Matrixx::getWidth() const
	{
		return width;
	}

	Roww Matrixx::getRow(const int row) const
	{
		Roww copyRow(row);
		for (int col = 0; col < width; col++) {
			copyRow[col] = rows[row][col];
		}
		return copyRow;
	}
	Vectorr Matrixx::getColumn(const int col) const
	{
		Vectorr copyVector(width);
		for (int row = 0; row < height; row++) {
			copyVector[row] = rows[row][col];
		}
		return copyVector;
	}

	const std::string Matrixx::str() const
	{
		std::string matrixString = "(" + std::to_string(height) + " x " + std::to_string(width) + " Matrix)\n";
		for (int row = 0; row < height; row++) {
			matrixString += rows[row].str();
		}
		return matrixString;
	}






	Roww::Roww(const int width)
	{
		init(width);
	}
	Roww::Roww(const Roww& copyRow)
		: Roww(copyRow.width)
	{
		for (int col = 0; col < width; col++) {
			entries[col] = copyRow[col];
		}
	}
	linalg::Roww::~Roww()
	{
		delete[] entries;
	}
	void Roww::init(int width)
	{
		check::LengthInfo lengthInfo{ 1, width };
		int errorNumber = check::checkWidth(width);
		check::handleLengthError(errorNumber, lengthInfo);

		this->width = width;
		entries = new double[width];
	}

	double& Roww::operator[](const int col)
	{
		check::LengthInfo lengthInfo{ 1, width };
		check::IndexInfo indexInfo{ 0, col, lengthInfo };
		int errorNumber = check::checkColumnIndex(col, width);
		check::handleOutOfRange(errorNumber, indexInfo);

		return entries[col];
	}
	const double& Roww::operator[](const int col) const
	{
		check::LengthInfo lengthInfo{ 1, width };
		check::IndexInfo indexInfo{ 0, col, lengthInfo };
		int errorNumber = check::checkColumnIndex(col, width);
		check::handleOutOfRange(errorNumber, indexInfo);

		return entries[col];
	}

	Allocator& Roww::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Roww::allocate(const int sequence, const double value)
	{
		if (sequence < width) {
			entries[sequence] = check::preventNegativeZero(value);
		}
	}

	Roww Roww::operator+() const
	{
		return Roww(*this);
	}
	Roww Roww::operator-() const
	{
		Roww negativeRow(width);
		for (int col = 0; col < width; col++) {
			negativeRow[col] = check::preventNegativeZero(-entries[col]);
		}
		return negativeRow;
	}

	Roww& Roww::operator=(const Roww& rightRow)
	{
		if (this == &rightRow) {
			return *this;
		}

		Roww copyRow(rightRow);
		swap(*this, copyRow);
		return *this;
	}
	void swap(Roww& leftRow, Roww& rightRow) noexcept
	{
		std::swap(leftRow.entries, rightRow.entries);
		std::swap(leftRow.width, rightRow.width);
	}

	Roww& Roww::operator+=(const Roww& rightRow)
	{
		check::LengthInfo leftLengthInfo{ 1, width };
		check::LengthInfo rightLengthInfo{ 1, rightRow.width };
		check::OperationInfo operationInfo{ '+', leftLengthInfo, rightLengthInfo };
		int errorNumber = check::checkWidth(width, rightRow.width);
		check::handleLogicError(errorNumber, operationInfo);
		
		for (int col = 0; col < width; col++) {
			entries[col] = check::preventNegativeZero(entries[col] + rightRow[col]);
		}
		return *this;
	}
	Roww& Roww::operator-=(const Roww& rightRow)
	{
		check::LengthInfo leftLengthInfo{ 1, width };
		check::LengthInfo rightLengthInfo{ 1, rightRow.width };
		check::OperationInfo operationInfo{ '-', leftLengthInfo, rightLengthInfo };
		int errorNumber = check::checkWidth(width, rightRow.width);
		check::handleLogicError(errorNumber, operationInfo);
		
		for (int col = 0; col < width; col++) {
			entries[col] = check::preventNegativeZero(entries[col] - rightRow[col]);
		}
		return *this;
	}
	Roww& Roww::operator*=(const double multiplier)
	{
		for (int col = 0; col < width; col++) {
			entries[col] = check::preventNegativeZero(multiplier * entries[col]);
		}
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
	std::ostream& operator<<(std::ostream& outputStream, const Roww& outputRow)
	{
		outputStream << outputRow.str();
		return outputStream;
	}

	const int Roww::getWidth() const
	{
		return width;
	}

	const std::string Roww::str() const
	{
		std::string rowString = "[\t";
		for (int col = 0; col < width; col++) {
			rowString += std::to_string(entries[col]) + "\t";
		}
		rowString += "]\n";
		return rowString;
	}

	
	


	linalg::Vectorr::Vectorr(const int height)
	{
		init(height);
	}
	Vectorr::Vectorr(const Vectorr& copyVector)
		: Vectorr(copyVector.height)
	{
		for (int row = 0; row < height; row++) {
			entries[row] = copyVector[row];
		}
	}
	linalg::Vectorr::~Vectorr()
	{
		delete[] entries;
	}
	void Vectorr::init(const int height)
	{
		check::LengthInfo lengthInfo{ height, 1 };
		int errorNumber = check::checkHeight(height);
		check::handleLengthError(errorNumber, lengthInfo);

		this->height = height;
		entries = new double[height];
	}

	double& Vectorr::operator[](const int row)
	{
		check::LengthInfo lengthInfo{ height, 1 };
		check::IndexInfo indexInfo{ row, 0, lengthInfo };
		int errorNumber = check::checkRowIndex(row, height);
		check::handleOutOfRange(errorNumber, indexInfo);

		return entries[row];
	}
	const double& Vectorr::operator[](const int row) const
	{
		check::LengthInfo lengthInfo{ height, 1 };
		check::IndexInfo indexInfo{ row, 0, lengthInfo };
		int errorNumber = check::checkRowIndex(row, height);
		check::handleOutOfRange(errorNumber, indexInfo);

		return entries[row];
	}

	Allocator& Vectorr::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Vectorr::allocate(const int sequence, const double value)
	{
		if (sequence < height) {
			entries[sequence] = check::preventNegativeZero(value);
		}
	}

	Vectorr Vectorr::operator+() const
	{
		return Vectorr(*this);
	}
	Vectorr Vectorr::operator-() const
	{
		Vectorr negativeVector(height);
		for (int row = 0; row < height; row++) {
			negativeVector[row] = check::preventNegativeZero(-entries[row]);
		}
		return negativeVector;
	}

	Vectorr& Vectorr::operator=(const Vectorr& rightVector)
	{
		if (this == &rightVector) {
			return *this;
		}

		Vectorr copyVector(rightVector);
		swap(*this, copyVector);
		return *this;
	}
	void swap(Vectorr& leftVector, Vectorr& rightVector) noexcept
	{
		std::swap(leftVector.entries, rightVector.entries);
		std::swap(leftVector.height, rightVector.height);
	}

	Vectorr& Vectorr::operator+=(const Vectorr& rightVector)
	{
		check::LengthInfo leftLengthInfo{ height, 1 };
		check::LengthInfo rightLengthInfo{ rightVector.height, 1 };
		check::OperationInfo operationInfo{ '+', leftLengthInfo, rightLengthInfo };
		int errorNumber = check::checkHeight(height, rightVector.height);
		check::handleLogicError(errorNumber, operationInfo);
		
		for (int row = 0; row < height; row++) {
			entries[row] = check::preventNegativeZero(entries[row] + rightVector[row]);
		}
		return *this;
	}
	Vectorr& Vectorr::operator-=(const Vectorr& rightVector)
	{
		check::LengthInfo leftLengthInfo{ height, 1 };
		check::LengthInfo rightLengthInfo{ rightVector.height, 1 };
		check::OperationInfo operationInfo{ '-', leftLengthInfo, rightLengthInfo };
		int errorNumber = check::checkHeight(height, rightVector.height);
		check::handleLogicError(errorNumber, operationInfo);
		
		for (int row = 0; row < height; row++) {
			entries[row] = check::preventNegativeZero(entries[row] - rightVector[row]);
		}
		return *this;
	}
	Vectorr& Vectorr::operator*=(const double multiplier)
	{
		for (int row = 0; row < height; row++) {
			entries[row] = check::preventNegativeZero(multiplier * entries[row]);
		}
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
	std::ostream& operator<<(std::ostream& outputStream, const Vectorr& outputVector)
	{
		outputStream << outputVector.str();
		return outputStream;
	}

	const int Vectorr::getHeight() const
	{
		return height;
	}

	const std::string Vectorr::str() const
	{
		std::string vectorString = "(" + std::to_string(height) + " row Vector)\n";
		for (int row = 0; row < height; row++) {
			vectorString += "[\t" + std::to_string(entries[row]) + "\t]\n";
		}
		return vectorString;
	}

	


	Matrixx identityMatrix(const int length)
	{
		Matrixx identityMatrix(length, length);
		for (int row = 0; row < length; row++) {
			for (int col = 0; col < length; col++) {
				identityMatrix[row][col] = (row == col) ? 1.0 : 0.0;
			}
		}
		return identityMatrix;
	}

	Matrixx zeroMatrix(const int height, const int width)
	{
		Matrixx identityMatrix(height, width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				identityMatrix[row][col] = 0.0;
			}
		}
		return identityMatrix;
	}
	
	
	
	
	namespace check {
		const int checkHeight(const int height)
		{
			if (height < 1) {
				return (int)LengthState::InvalidHeight;
			}
			return (int)LengthState::NoExcept;
		}
		const int checkWidth(const int width)
		{
			if (width < 1) {
				return (int)LengthState::InvalidWidth;
			}
			return (int)LengthState::NoExcept;
		}
		void handleLengthError(const int errorNumber, const LengthInfo lengthInfo)
		{
			switch (errorNumber) {
			case (int)LengthState::InvalidHeight:
				throw std::length_error("Bad height : " + std::to_string(lengthInfo.height) + "\n");
			case (int)LengthState::InvalidWidth:
				throw std::length_error("Bad width : " + std::to_string(lengthInfo.width) + "\n");
			case (int)LengthState::InvalidHeightAndWidth:
				throw std::length_error(
					"Bad height : " + std::to_string(lengthInfo.height) + "\n"
					+ "Bad width : " + std::to_string(lengthInfo.width) + "\n");
			default:
				return;
			}
		}

		const int checkRowIndex(const int row, const int height)
		{
			if (row >= height || row < 0) {
				return (int)IndexState::RowIndexOutOfRange;
			}
			return (int)IndexState::NoExcept;
		}
		const int checkColumnIndex(const int col, const int width)
		{
			if (col >= width || col < 0) {
				return (int)IndexState::ColumnIndexOutOfRange;
			}
			return (int)IndexState::NoExcept;
		}
		void handleOutOfRange(const int errorNumber, IndexInfo indexInfo)
		{
			switch (errorNumber) {
			case (int)IndexState::RowIndexOutOfRange:
				throw std::out_of_range(
					"Row index out of range(0-" + std::to_string(indexInfo.lengthInfo.height - 1)
					+ ") : " + std::to_string(indexInfo.row) + "\n");
			case (int)IndexState::ColumnIndexOutOfRange:
				throw std::out_of_range(
					"Column index out of range(0-" + std::to_string(indexInfo.lengthInfo.width - 1)
					+ ") : " + std::to_string(indexInfo.col) + "\n");
			case (int)IndexState::BothIndexOutOfRange:
				throw std::out_of_range(
					"Row index out of range(0-" + std::to_string(indexInfo.lengthInfo.height - 1)
					+ ") : " + std::to_string(indexInfo.row) + "\n"
					+ "Column index out of range(0-" + std::to_string(indexInfo.lengthInfo.width - 1)
					+ ") : " + std::to_string(indexInfo.col) + "\n");
			default:
				return;
			}
		}

		const int checkHeight(const int height1, const int height2)
		{
			if (height1 != height2) {
				return (int)OperationState::HeightDoNotMatch;
			}
			return (int)OperationState::NoExcept;
		}
		const int checkWidth(const int width1, const int width2)
		{
			if (width1 != width2) {
				return (int)OperationState::WidthDoNotMatch;
			}
			return (int)OperationState::NoExcept;
		}
		const int checkJoinLength(const int width, const int height)
		{
			if (width != height) {
				return (int)OperationState::JoinLengthDoNotMatch;
			}
			return (int)OperationState::NoExcept;
		}
		void handleLogicError(const int errorNumber, OperationInfo operationInfo)
		{
			switch (errorNumber) {
			case (int)OperationState::HeightDoNotMatch:
				throw std::length_error(
					"Height do not match : ("
					+ std::to_string(operationInfo.lengthInfo1.height) + " x "
					+ std::to_string(operationInfo.lengthInfo1.width) + ") "
					+ operationInfo.operation + " ("
					+ std::to_string(operationInfo.lengthInfo2.height) + " x "
					+ std::to_string(operationInfo.lengthInfo2.width) + ")\n");
			case (int)OperationState::WidthDoNotMatch:
				throw std::length_error(
					"Width do not match : ("
					+ std::to_string(operationInfo.lengthInfo1.height) + " x "
					+ std::to_string(operationInfo.lengthInfo1.width) + ") "
					+ operationInfo.operation + " ("
					+ std::to_string(operationInfo.lengthInfo2.height) + " x "
					+ std::to_string(operationInfo.lengthInfo2.width) + ")\n");
			case (int)OperationState::BothLengthDoNotMatch:
				throw std::length_error(
					"Height and width do not match : ("
					+ std::to_string(operationInfo.lengthInfo1.height) + " x "
					+ std::to_string(operationInfo.lengthInfo1.width) + ") "
					+ operationInfo.operation + " ("
					+ std::to_string(operationInfo.lengthInfo2.height) + " x "
					+ std::to_string(operationInfo.lengthInfo2.width) + ")\n");
			case (int)OperationState::JoinLengthDoNotMatch:
				throw std::length_error(
					"Cannot multiply : ("
					+ std::to_string(operationInfo.lengthInfo1.height) + " x "
					+ std::to_string(operationInfo.lengthInfo1.width)
					+ ") * ("
					+ std::to_string(operationInfo.lengthInfo2.height) + " x "
					+ std::to_string(operationInfo.lengthInfo2.width) + ")\n");
			default:
				return;
			}
		}

		double preventNegativeZero(const double value)
		{
			return (value == -0.0) ? 0.0 : value;
		}


	}
}