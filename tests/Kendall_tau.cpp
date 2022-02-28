// C++ program for Merge Sort
#include <iostream>
#include <vector>
using namespace std;
  
// // Merges two subarrays of array[].
// // First subarray is arr[begin..mid]
// // Second subarray is arr[mid+1..end]
// void merge(int array[], int const left, int const mid, int const right)
// {
//     auto const subArrayOne = mid - left + 1;
//     auto const subArrayTwo = right - mid;
  
//     // Create temp arrays
//     auto *leftArray = new int[subArrayOne],
//          *rightArray = new int[subArrayTwo];
  
//     // Copy data to temp arrays leftArray[] and rightArray[]
//     for (auto i = 0; i < subArrayOne; i++)
//         leftArray[i] = array[left + i];
//     for (auto j = 0; j < subArrayTwo; j++)
//         rightArray[j] = array[mid + 1 + j];
  
//     auto indexOfSubArrayOne = 0, // Initial index of first sub-array
//         indexOfSubArrayTwo = 0; // Initial index of second sub-array
//     int indexOfMergedArray = left; // Initial index of merged array
  
//     // Merge the temp arrays back into array[left..right]
//     while (indexOfSubArrayOne < subArrayOne && indexOfSubArrayTwo < subArrayTwo) {
//         if (leftArray[indexOfSubArrayOne] <= rightArray[indexOfSubArrayTwo]) {
//             array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
//             indexOfSubArrayOne++;
//         }
//         else {
//             array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
//             indexOfSubArrayTwo++;
//         }
//         indexOfMergedArray++;
//     }
//     // Copy the remaining elements of
//     // left[], if there are any
//     while (indexOfSubArrayOne < subArrayOne) {
//         array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
//         indexOfSubArrayOne++;
//         indexOfMergedArray++;
//     }
//     // Copy the remaining elements of
//     // right[], if there are any
//     while (indexOfSubArrayTwo < subArrayTwo) {
//         array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
//         indexOfSubArrayTwo++;
//         indexOfMergedArray++;
//     }
// }
  
// // begin is for left index and end is
// // right index of the sub-array
// // of arr to be sorted */
// void mergeSort(int array[], int const begin, int const end)
// {
//     if (begin >= end)
//         return; // Returns recursively
  
//     auto mid = begin + (end - begin) / 2;
//     mergeSort(array, begin, mid);
//     mergeSort(array, mid + 1, end);
//     merge(array, begin, mid, end);
// }
  
// // UTILITY FUNCTIONS
// // Function to print an array
// void printArray(int A[], int size)
// {
//     for (auto i = 0; i < size; i++)
//         cout << A[i] << " ";
// }

// Merges two subarrays of array[].
// First subarray is arr[begin..mid]
// Second subarray is arr[mid+1..end]
size_t Kendall_tau_merge(int array[], int const left, int const mid, int const right)
{
    size_t num_inv = 0;
    auto const subArrayOne = mid - left + 1;
    auto const subArrayTwo = right - mid;
  
    // Create temp arrays
    auto *leftArray = new int[subArrayOne],
         *rightArray = new int[subArrayTwo];
  
    // Copy data to temp arrays leftArray[] and rightArray[]
    for (auto i = 0; i < subArrayOne; i++)
        leftArray[i] = array[left + i];
    for (auto j = 0; j < subArrayTwo; j++)
        rightArray[j] = array[mid + 1 + j];
  
    auto indexOfSubArrayOne = 0, // Initial index of first sub-array
        indexOfSubArrayTwo = 0; // Initial index of second sub-array
    int indexOfMergedArray = left; // Initial index of merged array
  
    // Merge the temp arrays back into array[left..right]
    while (indexOfSubArrayOne < subArrayOne && indexOfSubArrayTwo < subArrayTwo) {
        if (leftArray[indexOfSubArrayOne] <= rightArray[indexOfSubArrayTwo]) {
            array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
            indexOfSubArrayOne++;
        }
        else {
            array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
            num_inv += subArrayOne - indexOfSubArrayOne;
            indexOfSubArrayTwo++;
        }
        indexOfMergedArray++;
    }
    // Copy the remaining elements of
    // left[], if there are any
    while (indexOfSubArrayOne < subArrayOne) {
        array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
        indexOfSubArrayOne++;
        indexOfMergedArray++;
    }
    // Copy the remaining elements of
    // right[], if there are any
    while (indexOfSubArrayTwo < subArrayTwo) {
        array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
        indexOfSubArrayTwo++;
        indexOfMergedArray++;
    }
    return num_inv;
}
  
// Using merge sort to compute Kendal tau distance 
// array will be sorted in place
size_t Kendall_tau_inplace(int array[], int const begin, int const end)
{
    size_t num_inv = 0;
    if (begin >= end)
        return 0; // Returns recursively
    
    auto mid = begin + (end - begin) / 2;
    num_inv += Kendall_tau_inplace(array, begin, mid);
    num_inv += Kendall_tau_inplace(array, mid + 1, end);
    num_inv += Kendall_tau_merge(array, begin, mid, end);
    return num_inv;
}

template<typename T>
size_t Kendall_tau(const T& perm)
{
    // Create temp array
    size_t len = perm.size();
    auto *arr = new int[len];
    // Copy data to temp array
    for (auto i = 0; i < len; i++)
        arr[i] = perm[i];
    size_t num_inv = Kendall_tau_inplace(arr, 0, len - 1);

    return num_inv;
}


// Driver code
int main()
{
    int arr[] = { 12, 11, 13, 5, 6, 7 };
    auto arr_size = sizeof(arr) / sizeof(arr[0]);
  
    // cout << "Given array is \n";
    // printArray(arr, arr_size);
    // std::cout << "\n";
    // // mergeSort(arr, 0, arr_size - 1);
    // std::cout << "Kendall tau distance is " << Kendall_tau(arr) << std::endl;
    // cout << "\nSorted array is \n";
    // printArray(arr, arr_size);
    // std::cout << "\n";

    std::vector<int> vec (arr, arr + arr_size ) ;
    std::cout << "Kendall tau distance is " << Kendall_tau(vec) << std::endl;
    return 0;
}