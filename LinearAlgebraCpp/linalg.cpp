#include "linalg.h"

namespace linalg {

	Matrixx::Matrixx(const int height, const int width)
		: Tensorr(), impl(std::make_unique<Impl>(height, width))
	{
	}
	Matrixx::Matrixx(const Matrixx& copyMatrix)
		: Tensorr(), impl(std::make_unique<Impl>(*(copyMatrix.impl)))
	{
	}
	/*Matrixx::Matrixx(Matrixx&& moveMatrix) noexcept
	{
		swap(*this, moveMatrix);
	}*/
	Matrixx::Matrixx(const Roww& copyRow)
		: Tensorr(), impl(std::make_unique<Impl>(*(copyRow.impl)))
	{
	}
	Matrixx::Matrixx(const Vectorr& copyVector)
		: Tensorr(), impl(std::make_unique<Impl>(*(copyVector.impl)))
	{
	}
	Matrixx::Matrixx(const Impl& matrixImpl)
		: Tensorr(), impl(std::make_unique<Impl>(matrixImpl))
	{
	}
	void Matrixx::init(const int height, const int width)
	{		
		impl->init(height, width);
	}

	void Matrixx::reduce()
	{
		impl->reduce();
	}

	void Matrixx::toEchelonForm()
	{
		impl->toEchelonForm();
	}

	void Matrixx::toReducedEchelonForm()
	{
		impl->toReducedEchelonForm();
	}

	bool Matrixx::isEchelonForm()
	{
		return impl->isEchelonForm();
	}

	Matrixx Matrixx::block(const size_t beginRow, const size_t beginCol,
		const size_t blockHeight, const size_t blockWidth) const
	{
		return impl->block(beginRow, beginCol, blockHeight, blockWidth);
	}

	Matrixx Matrixx::inverse()
	{
		return impl->inverse();
	}

	Matrixx Matrixx::transpose(const bool inplace)
	{
		Matrixx transposedMatrix = impl->transpose();
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
		return Impl::identity(length);
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

	const Roww& Matrixx::operator[](const size_t row) const
	{
		return (*impl)[row];
	}
	Roww& Matrixx::operator[](const size_t row)
	{
		return const_cast<Roww&>(static_cast<const Matrixx&>(*this)[row]);
	}

	const Roww& Matrixx::operator()(const int row) const
	{
		return (*impl)(row);
	}
	Roww& Matrixx::operator()(const int row)
	{
		return const_cast<Roww&>(static_cast<const Matrixx&>(*this)(row));
	}

	const double& Matrixx::operator()(const int row, const int col) const
	{
		return (*impl)(row, col);
	}
	double& Matrixx::operator()(const int row, const int col)
	{
		return const_cast<double&>(static_cast<const Matrixx&>(*this)(row, col));
	}

	
	Allocatorr& Matrixx::operator<<(const double value)
	{
		impl->allocate(0, value);
		return *(new Allocatorr(*this, 1));
	}
	void Matrixx::allocate(const size_t sequence, const double value)
	{
		impl->allocate(sequence, value);
	}
	Matrixx& Matrixx::operator=(const std::initializer_list<double> values)
	{
		impl->allocate(values);
		return *this;
	}

	Matrixx Matrixx::operator+() const
	{
		return +(*impl);
	}
	Matrixx Matrixx::operator-() const
	{
		return -(*impl);
	}

	void swap(Matrixx& leftMatrix, Matrixx& rightMatrix) noexcept
	{
		std::swap(leftMatrix.impl, rightMatrix.impl);
	}
	Matrixx& Matrixx::operator=(const Matrixx& rightMatrix)
	{
		if (this == &rightMatrix) {
			return *this;
		}

		*impl = *(rightMatrix.impl);
		return *this;
	}
	/*Matrixx& Matrixx::operator=(Matrixx&& rightMatrix)
	{
		Matrixx moveMatrix(std::move(rightMatrix));
		swap(*this, moveMatrix);
		return *this;
	}*/
	Matrixx& Matrixx::operator+=(const Matrixx& rightMatrix)
	{
		*impl += *(rightMatrix.impl);
		return *this;
	}
	Matrixx& Matrixx::operator-=(const Matrixx& rightMatrix)
	{
		* impl -= *(rightMatrix.impl);
		return *this;
	}
	Matrixx& Matrixx::operator*=(const double multiplier)
	{
		*impl *= multiplier;
		return *this;
	}
	Matrixx& Matrixx::operator*=(const Matrixx& rightMatrix)
	{
		*impl *= *(rightMatrix.impl);
		return *this;
	}
	Matrixx& Matrixx::operator/=(const double divisor)
	{
		*impl /= divisor;
		return *this;
	}

	Matrixx& Matrixx::operator&=(const Matrixx& rightMatrix)
	{
		*impl &= *(rightMatrix.impl);
		return *this;
	}
	Matrixx& Matrixx::operator&=(const Vectorr& rightVector)
	{
		*impl &= *(rightVector.impl);
		return *this;
	}
	Matrixx& Matrixx::operator|=(const Matrixx& lowerMatrix)
	{
		*impl |= *(lowerMatrix.impl);
		return *this;
	}
	Matrixx& Matrixx::operator|=(const Roww& lowerRow)
	{
		*impl |= *(lowerRow.impl);
		return *this;
	}

	Matrixx operator+(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		return *(leftMatrix.impl) + *(rightMatrix.impl);
	}
	Matrixx operator-(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		return *(leftMatrix.impl) - *(rightMatrix.impl);
	}
	Matrixx operator*(const double multiplier, const Matrixx& rightMatrix)
	{
		return *(rightMatrix.impl) * multiplier;
	}
	Matrixx operator*(const Matrixx& leftMatrix, const double multiplier)
	{
		return multiplier * leftMatrix;
	}
	Matrixx operator*(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		return *(leftMatrix.impl) * *(rightMatrix.impl);
	}
	Vectorr operator*(const Matrixx& leftMatrix, const Vectorr& rightVector)
	{
		return Vectorr(*(leftMatrix.impl) * *(rightVector.impl));
	}
	Matrixx operator/(const Matrixx& leftMatrix, const double divisor)
	{
		return *(leftMatrix.impl) / divisor;
	}

	Matrixx operator&(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		return *(leftMatrix.impl) & *(rightMatrix.impl);
	}
	Matrixx operator&(const Matrixx& leftMatrix, const Vectorr& rightVector)
	{
		return *(leftMatrix.impl) & *(rightVector.impl);
	}
	Matrixx operator&(const Vectorr& leftVector, const Matrixx& rightMatrix)
	{
		return *(leftVector.impl) & *(rightMatrix.impl);
	}
	Matrixx operator&(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		return *(leftVector.impl) & *(rightVector.impl);
	}

	Matrixx operator|(const Matrixx& upperMatrix, const Matrixx& lowerMatrix)
	{
		return *(upperMatrix.impl) | *(lowerMatrix.impl);
	}
	Matrixx operator|(const Matrixx& upperMatrix, const Roww& lowerRow)
	{
		return *(upperMatrix.impl) | *(lowerRow.impl);
	}
	Matrixx operator|(const Roww& upperRow, const Matrixx& lowerMatrix)
	{
		return *(upperRow.impl) | *(lowerMatrix.impl);
	}
	Matrixx operator|(const Roww& upperRow, const Roww& lowerRow)
	{
		return *(upperRow.impl) | *(lowerRow.impl);
	}

	bool operator==(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		return *(leftMatrix.impl) == *(rightMatrix.impl);
	}
	bool operator!=(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		return *(leftMatrix.impl) != *(rightMatrix.impl);
	}

	std::ostream& operator<<(std::ostream& outputStream, const Matrixx& outputMatrix)
	{
		outputStream << outputMatrix.str();
		return outputStream;
	}

	const size_t Matrixx::size() const
	{
		return impl->height() * impl->width();
	}

	const size_t Matrixx::height() const
	{
		return impl->height();
	}
	const size_t Matrixx::width() const
	{
		return impl->width();
	}

	Roww Matrixx::getRow(const int row) const
	{
		return impl->getRow(row);
	}
	Vectorr Matrixx::getColumn(const int col) const
	{
		return impl->getColumn(col);
	}

	const std::string Matrixx::str() const
	{
		return impl->str();
	}





	Roww::Roww(const int size)
		: Tensorr(), impl(std::make_unique<Impl>(size))
	{
	}
	Roww::Roww(const Roww& copyRow)
		: Tensorr(), impl(std::make_unique<Impl>(*(copyRow.impl)))
	{
	}
	Roww::Roww(const Impl& rowImpl)
		: Tensorr(), impl(std::make_unique<Impl>(rowImpl))
	{
	}
	/*Roww::Roww(Roww&& moveRow) noexcept
	{
		swap(*this, moveRow);
	}*/
	void Roww::init(const int size)
	{
		impl->init(size);
	}

	const double& Roww::operator[](const size_t col) const
	{
		return (*impl)[col];
	}
	double& Roww::operator[](const size_t col)
	{
		return const_cast<double&>(static_cast<const Roww&>(*this)[col]);
	}

	const double& Roww::operator()(const int col) const
	{
		return (*impl)(col);
	}
	double& Roww::operator()(const int col)
	{
		return const_cast<double&>(static_cast<const Roww&>(*this)(col));
	}

	Allocatorr& Roww::operator<<(const double value)
	{
		impl->allocate(0, value);
		return *(new Allocatorr(*this, 1));
	}
	void Roww::allocate(const size_t sequence, const double value)
	{
		impl->allocate(sequence, value);
	}
	Roww& Roww::operator=(const std::initializer_list<double> values)
	{
		impl->allocate(values);
		return *this;
	}

	Roww Roww::operator+() const
	{
		return +(*impl);
	}
	Roww Roww::operator-() const
	{
		return -(*impl);
	}

	void swap(Roww& leftRow, Roww& rightRow) noexcept
	{
		std::swap(leftRow.impl, rightRow.impl);
	}
	Roww& Roww::operator=(const Roww& rightRow)
	{
		if (this == &rightRow) {
			return *this;
		}

		*impl = *(rightRow.impl);
		return *this;
	}
	/*Roww& Roww::operator=(Roww&& rightRow) noexcept
	{
		Roww moveRow(std::move(rightRow));
		swap(*this, rightRow);
		return *this;
	}*/
	Roww& Roww::operator+=(const Roww& rightRow)
	{
		*impl += *(rightRow.impl);
		return *this;
	}
	Roww& Roww::operator-=(const Roww& rightRow)
	{
		*impl -= *(rightRow.impl);
		return *this;
	}
	Roww& Roww::operator*=(const double multiplier)
	{
		*impl *= multiplier;
		return *this;
	}
	Roww& Roww::operator/=(const double divisor)
	{
		*impl /= divisor;
		return *this;
	}

	Roww& Roww::operator&=(const Roww& rightRow)
	{
		*impl &= *(rightRow.impl);
		return *this;
	}

	Roww operator+(const Roww& leftRow, const Roww& rightRow)
	{
		return *(leftRow.impl) + *(rightRow.impl);
	}
	Roww operator-(const Roww& leftRow, const Roww& rightRow)
	{
		return *(leftRow.impl) - *(rightRow.impl);
	}
	Roww operator*(const double multiplier, const Roww& rightRow)
	{
		return *(rightRow.impl) * multiplier;
	}
	Roww operator*(const Roww& leftRow, const double multiplier)
	{
		return *(leftRow.impl) * multiplier;
	}
	Roww operator/(const Roww& leftRow, const double divisor)
	{
		return *(leftRow.impl) / divisor;
	}
	Roww operator&(const Roww& leftRow, const Roww& rightRow)
	{
		return *(leftRow.impl) & *(rightRow.impl);
	}
	bool operator==(const Roww& leftRow, const Roww& rightRow)
	{
		return *(leftRow.impl) == *(rightRow.impl);
	}
	bool operator!=(const Roww& leftRow, const Roww& rightRow)
	{
		return *(leftRow.impl) != *(rightRow.impl);
	}
	std::ostream& operator<<(std::ostream& outputStream, const Roww& outputRow)
	{
		outputStream << outputRow.str();
		return outputStream;
	}

	const size_t Roww::size() const
	{
		return impl->width();
	}
	const size_t Roww::width() const
	{
		return impl->width();
	}

	const std::string Roww::str() const
	{
		return impl->str();
	}





	Vectorr::Vectorr(const int size)
		: Tensorr(), impl(std::make_unique<Impl>(size))
	{
	}
	Vectorr::Vectorr(const Vectorr& copyVector)
		: Tensorr(), impl(std::make_unique<Impl>(*(copyVector.impl)))
	{
	}
	Vectorr::Vectorr(const Impl& vectorImpl)
		: Tensorr(), impl(std::make_unique<Impl>(vectorImpl))
	{
	}
	/*Vectorr::Vectorr(Vectorr&& moveVector) noexcept
	{
		swap(*this, moveVector);
	}*/
	void Vectorr::init(const int size)
	{
		impl->init(size);
	}

	const double& Vectorr::operator[](const size_t row) const
	{
		return (*impl)[row];
	}
	double& Vectorr::operator[](const size_t row)
	{
		return const_cast<double&>(static_cast<const Vectorr&>(*this)[row]);
	}

	const double& Vectorr::operator()(const int row) const
	{
		return (*impl)(row);
	}
	double& Vectorr::operator()(const int row)
	{
		return const_cast<double&>(static_cast<const Vectorr&>(*this)(row));
	}

	Allocatorr& Vectorr::operator<<(const double value)
	{
		impl->allocate(0, value);
		return *(new Allocatorr(*this, 1));
	}
	void Vectorr::allocate(const size_t sequence, const double value)
	{
		impl->allocate(sequence, value);
	}
	Vectorr& Vectorr::operator=(const std::initializer_list<double> values)
	{
		impl->allocate(values);
		return *this;
	}
	

	Vectorr Vectorr::operator+() const
	{
		return +(*impl);
	}
	Vectorr Vectorr::operator-() const
	{
		return -(*impl);
	}

	void swap(Vectorr& leftVector, Vectorr& rightVector) noexcept
	{
		std::swap(leftVector.impl, rightVector.impl);
	}
	Vectorr& Vectorr::operator=(const Vectorr& rightVector)
	{
		if (this == &rightVector) {
			return *this;
		}

		*impl = *(rightVector.impl);
		return *this;
	}
	/*Vectorr& Vectorr::operator=(Vectorr&& rightVector) noexcept
	{
		Vectorr moveVector(std::move(rightVector));
		swap(*this, moveVector);
		return *this;
	}*/
	Vectorr& Vectorr::operator+=(const Vectorr& rightVector)
	{
		*impl += *(rightVector.impl);
		return *this;
	}
	Vectorr& Vectorr::operator-=(const Vectorr& rightVector)
	{
		*impl -= *(rightVector.impl);
		return *this;
	}
	Vectorr& Vectorr::operator*=(const double multiplier)
	{
		*impl *= multiplier;
		return *this;
	}
	Vectorr& Vectorr::operator/=(const double divisor)
	{
		*impl /= divisor;
		return *this;
	}

	Vectorr& Vectorr::operator|=(const Vectorr& lowerVector)
	{
		*impl |= *(lowerVector.impl);
		return *this;
	}

	Vectorr operator+(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		return *(leftVector.impl) + *(rightVector.impl);
	}
	Vectorr operator-(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		return *(leftVector.impl) - *(rightVector.impl);
	}
	Vectorr operator*(const double multiplier, const Vectorr& rightVector)
	{
		return *(rightVector.impl) * multiplier;
	}
	Vectorr operator*(const Vectorr& leftVector, const double multiplier)
	{
		return *(leftVector.impl) * multiplier;
	}
	
	Vectorr operator/(const Vectorr& leftVector, const double divisor)
	{
		return *(leftVector.impl) / divisor;
	}

	Vectorr operator|(const Vectorr& upperVector, const Vectorr& lowerVector)
	{
		return *(upperVector.impl) | *(lowerVector.impl);
	}
	bool operator==(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		return *(leftVector.impl) == *(rightVector.impl);
	}
	bool operator!=(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		return *(leftVector.impl) != *(rightVector.impl);
	}
	std::ostream& operator<<(std::ostream& outputStream, const Vectorr& outputVector)
	{
		outputStream << outputVector.str();
		return outputStream;
	}

	const size_t Vectorr::size() const
	{
		return impl->height();
	}
	const size_t Vectorr::height() const
	{
		return impl->height();
	}

	const std::string Vectorr::str() const
	{
		return impl->str();
	}
}