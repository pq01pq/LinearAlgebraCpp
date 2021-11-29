#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <memory>

namespace linalg {
	class Matrixx;
	class Roww;
	class Vectorr;
	//class Celll;
	class Allocatable;
	class Allocator;

	class Allocator {
		friend class Allocatable;
	public:
		Allocator(Allocatable& target, const int sequence);
		~Allocator();
		Allocator& operator,(const double value);
		
	private:
		Allocatable& target;
		int sequence;
	};

	class Allocatable {
		friend class Allocator;
	protected:
		virtual void allocate(const int sequence, const double value) {};
	};

	class Matrixx : private Allocatable {
		friend class Roww;
		friend class Vectorr;
		//friend class Celll;
	public:
		Matrixx(int height, int width);
		Matrixx(const Matrixx& copyMatrix);
		~Matrixx();

		Roww& operator[](int row);
		const Roww& operator[](int row) const;
		double& operator()(int row, int col) const;

		Allocator& operator<<(const double value);

		Matrixx& operator=(const Matrixx& rightMatrix);
		

		const int getHeight() const;
		const int getWidth() const;

		Roww& getRow(int row) const;
		Vectorr& getColumn(int col) const;

		const std::string str() const;

	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		Roww* rows;
		int height, width;

		struct Pivot {
			int row, col;
			double value;
		};

		friend void swap(Matrixx& leftMatrix, Matrixx& rightMatrix) noexcept;
	};

	class Roww : private Allocatable {
		friend class Matrixx;
		friend class Vectorr;
		//friend class Celll;
	public:
		Roww() = default;
		explicit Roww(int width);
		Roww(const Roww& copyRow);
		~Roww();

		double& operator[](int col);
		const double& operator[](int col) const;

		Allocator& operator<<(const double value);

		Roww& operator=(const Roww& rightRow);

		const int getWidth() const;

		const std::string str() const;
	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		double* entries;
		int width;

		void init(int width);

		friend void swap(Roww& leftRow, Roww& rightRow) noexcept;
	};

	class Vectorr : private Allocatable {
		friend class Matrixx;
		friend class Roww;
		//friend class Celll;
	public:
		Vectorr() = default;
		explicit Vectorr(int height);
		Vectorr(const Vectorr& copyVector);
		~Vectorr();

		double& operator[](int row);
		const double& operator[](int row) const;

		Allocator& operator<<(const double value);

		Vectorr& operator=(const Vectorr& rightVector);
		
		const int getHeight() const;

		const std::string str() const;
	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		double* entries;
		int height;

		void init(int height);

		friend void swap(Vectorr& leftVector, Vectorr& rightVector) noexcept;
	};

	//class Celll {
	//	friend class Matrixx;
	//	friend class Roww;
	//	friend class Vectorr;
	//public:
	//	Celll() = default;
	//	Celll(double value);
	//	//Celll(const Celll& copyCell);
	//	~Celll();

	//	operator double() const;
	//	//Celll& operator=(const Celll& rightCell);

	//	void set(double value);
	//	const double get() const;
	//private:
	//	double value;
	//};
}
