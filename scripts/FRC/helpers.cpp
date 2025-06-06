#include <assert.h>
#include "helpers.h"

void bounds_check(int len, int start, int end){
	assert(start >= 0 && end < len && start <= end);
}