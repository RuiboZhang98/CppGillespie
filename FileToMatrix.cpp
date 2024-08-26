/*
Copyright (c) 2020, Daniel Hanson
CFRM 520, University of Washington
For instructional purposes only
*/

#include "FileToMatrix.h"

using std::vector;
using std::fstream;
using std::ifstream;
using std::string;

using std::cout;		// for debug mode
using std::endl;		// for debug mode

// References:
//http://ubuntuforums.org/showthread.php?t=2002692
//http://www.cplusplus.com/reference/iostream/ifstream/
//http://www.cplusplus.com/reference/string/getline/

FileToMatrix::FileToMatrix(long nrow, long ncol, const std::string& filePathAndName) :nrow_(nrow), ncol_(ncol), 
	filePathAndName_(filePathAndName)
{
	csvToMatrix_();
}

Eigen::MatrixXd FileToMatrix::operator()() const
{
	return eigDataMatrix_;
}

void FileToMatrix::csvToMatrix_()
{
	eigDataVector_.resize(nrow_ * ncol_);
	dataAsString_.clear();

	ifstream file(filePathAndName_);
	double d;
	string s;

	while(file.good())
	{
		getline(file, s);

		std::stringstream lineStream(s);
        std::string cell;
        while(std::getline(lineStream, cell,','))
        {
			dataAsString_.push_back(cell);
        }
	}

	for (long k = 0; k < long(dataAsString_.size()); ++k)	// this is not elegant as using iterators, but we don't have this feature in Eigen
	{
		d = strtod(dataAsString_.at(k).c_str(), NULL);
		eigDataVector_(k) = d;
	}

	constructEigDataMatrix_();
}

void FileToMatrix::constructEigDataMatrix_()
{
	eigDataMatrix_ = eigDataVector_;
	eigDataMatrix_.resize(ncol_, nrow_);		// This is how I got Eigen to put the data in row-wise.  RowVector didn't seem to work.
	eigDataMatrix_.transposeInPlace();			// Now transpose the matrix to get the desired result.
}

/*
 Copyright 2020 <Daniel Hanson  MIT License

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
 subject to the following conditions:

 The above copyright notice and this permission
 notice shall be included in all copies
 or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
 FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
