#ifndef __QSORT__
#define __QSORT__


template<class TYPE, class SIZE> void QuickSort1(SIZE len, const TYPE *data, SIZE *order, int direction = 0);
template<class TYPE, class SIZE> void QuickSort2(SIZE len, TYPE *data, int direction = 0);


///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class TYPE, class SIZE> void QuickSortkernel1(const TYPE *data, SIZE *order, SIZE low, SIZE high, int direction);
template<class TYPE, class SIZE> void QuickSortkernel2(TYPE *data, SIZE low, SIZE high, int direction);

template<class TYPE, class SIZE> void QuickSort1(SIZE len, const TYPE *data, SIZE *order, int direction)
{
	for (SIZE i = 0; i < len; ++i)
	{
		order[i] = i;
	}
	QuickSortkernel1(data, order, (SIZE)0, len - 1, direction);
}

template<class TYPE, class SIZE> void QuickSort2(SIZE len, TYPE *data, int direction)
{
	QuickSortkernel2(data, (SIZE)0, len - 1, direction);
}

template<class TYPE, class SIZE> void QuickSortkernel1(const TYPE *data, SIZE *order, SIZE start, SIZE end, int direction)
{
	if (0 == direction)
	{
		while (start < end)
		{
			SIZE i = start;
			SIZE j = end + 1;
			SIZE mid = (i + j) >> 1;
			SIZE o = order[mid];
			order[mid] = order[start];
			order[start] = o;
			TYPE x = data[o];

			do
			{
				do ++i; while (i < end && data[order[i]] < x);
				do --j; while (j > start && data[order[j]] > x);
				if (i < j)
				{
					SIZE o2 = order[i];
					order[i] = order[j];
					order[j] = o2;
				}
			} while (i < j);

			order[start] = order[j];
			order[j] = o;

			if (j - start > end - j)
			{
				QuickSortkernel1(data, order, j + 1, end, direction);
				end = j - 1;
			}
			else
			{
				QuickSortkernel1(data, order, start, j - 1, direction);
				start = j + 1;
			}
		}
	}
	else
	{
		while (start < end)
		{
			SIZE i = start;
			SIZE j = end + 1;
			SIZE mid = (i + j) >> 1;
			SIZE o = order[mid];
			order[mid] = order[start];
			order[start] = o;
			TYPE x = data[o];

			do
			{
				do ++i; while (i < end && data[order[i]] > x);
				do --j; while (j > start && data[order[j]] < x);
				if (i < j)
				{
					SIZE o2 = order[i];
					order[i] = order[j];
					order[j] = o2;
				}
			} while (i < j);

			order[start] = order[j];
			order[j] = o;

			if (j - start > end - j)
			{
				QuickSortkernel1(data, order, j + 1, end, direction);
				end = j - 1;
			}
			else
			{
				QuickSortkernel1(data, order, start, j - 1, direction);
				start = j + 1;
			}
		}
	}
}

template<class TYPE, class SIZE> void QuickSortkernel2(TYPE *data, SIZE start, SIZE end, int direction)
{
	if (0 == direction)
	{
		while (start < end)
		{
			SIZE i = start;
			SIZE j = end + 1;
			SIZE mid = (i + j) >> 1;
			TYPE x = data[mid];
			data[mid] = data[start];
			data[start] = x;

			do
			{
				do ++i; while (i < end && data[i] < x);
				do --j; while (j > start && data[j] > x);
				if (i < j)
				{
					TYPE x2 = data[i];
					data[i] = data[j];
					data[j] = x2;
				}
			} while (i < j);

			data[start] = data[j];
			data[j] = x;

			if (j - start > end - j)
			{
				QuickSortkernel2(data, j + 1, end, direction);
				end = j - 1;
			}
			else
			{
				QuickSortkernel2(data, start, j - 1, direction);
				start = j + 1;
			}
		}
	}
	else
	{
		while (start < end)
		{
			SIZE i = start;
			SIZE j = end + 1;
			SIZE mid = (i + j) >> 1;
			TYPE x = data[mid];
			data[mid] = data[start];
			data[start] = x;

			do
			{
				do ++i; while (i < end && data[i] > x);
				do --j; while (j > start && data[j] < x);
				if (i < j)
				{
					TYPE x2 = data[i];
					data[i] = data[j];
					data[j] = x2;
				}
			} while (i < j);

			data[start] = data[j];
			data[j] = x;

			if (j - start > end - j)
			{
				QuickSortkernel2(data, j + 1, end, direction);
				end = j - 1;
			}
			else
			{
				QuickSortkernel2(data, start, j - 1, direction);
				start = j + 1;
			}
		}
	}
}

#endif
