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
		
		const int getHeight() const;
		const int getWidth() const;

		const std::string str() const;
	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		Roww* rows;
		int height, width;
		int count = 0;

		struct Pivot {
			int row, col;
			double value;
		};
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

		const int getWidth() const;

		const std::string str() const;
	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		double* cells;
		int width;

		void init(int width);
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
		
		const int getHeight() const;

		const std::string str() const;
	protected:
		virtual void allocate(const int sequence, const double value) override;
	private:
		double* cells;
		int height;

		void init(int height);
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
