#include<bits/stdc++.h>
#include"E:\Matlab\extern\include\mex.h"

mxArray* ColumnTransform(mxArray* order) {
	const mwSize* mSize = mxGetDimensions(order);
	
	mxArray* ptr = (mxArray*)mxCalloc(*mSize, *mSize);

	if (*mSize < 1) return 0;
	for (int i = 0; i < *mSize; ++i) {
		*ptr(*(order + i) * (*mSize) + i) = 1;
	}
	return ((mxArray*)0);
}

