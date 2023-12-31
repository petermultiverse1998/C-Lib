//
// Created by peter on 6/21/2023.
//

#include "task.h"

#include <malloc.h>
#include <stdio.h>

static int allocatedMemory = 0;

typedef struct TaskQueueData TaskQueueData;

/**
 * It allocates the memory and return pointer to it
 * @param heap          : Pointer to static heap
 *                      : NULL for dynamic heap
 * @param sizeInByte    : Size in bytes
 * @return              : Pointer to allocated memory
 *                      : NULL if there exist no memory for allocation
 */
static void *allocateMemory(BuddyHeap *heap,int sizeInByte) {
    if(sizeInByte<=0)
        return NULL;
    void *ptr;
    ptr = heap != NULL ? StaticBuddyHeap.malloc(heap, sizeInByte) : malloc(sizeInByte);
    if (ptr != NULL)
        allocatedMemory += sizeInByte;
    return ptr;
}

/**
 * It free the allocated memory
 * @param heap          : Pointer to static heap
 *                      : NULL for dynamic heap
 * @param pointer       : Pointer to allocated Memory
 * @param sizeInByte    : Size to be freed
 * @return              : 1 for success (OR) 0 for failed
 */
static int freeMemory(BuddyHeap *heap,void *pointer, int sizeInByte) {
    if(pointer==NULL || sizeInByte<=0)
        return 0;
    heap != NULL ? StaticBuddyHeap.free(heap, pointer) : free(pointer);
    allocatedMemory -= sizeInByte;
    return 1;
}

/**
 * Check if memory pointer exist in defined memory
 */
static int validMemory(const char *func, BuddyHeap *heap, void *ptr) {
    if (ptr == NULL)
        return 1;
    uint64_t addr = (uint64_t) ptr;
    if (heap == NULL) {

        return 1;

        if (addr < 0x20000000 || addr > 0x20004fff) {
            printf("HashMap-%s:\n", func);
            printf("Memory : 0x20000000 - 0x20004fff\n");
            printf("Ptr : %p\n\n", ptr);
//			*(uint8_t*)NULL = 10;
            return 0;
        }
    } else {
        if (!StaticBuddyHeap.isValidPointer(*heap, ptr)) {
            printf("HashMap-%s:\n", func);
            printf("Heap : %p\n", heap);
            printf("Memory : %p - %p\n", heap->memory, heap->memory + heap->maxSize);
            printf("Ptr : %p\n\n", ptr);
//			*(uint8_t*) NULL = 10;
            return 0;
        }

    }
    return 1;
}

/**
 * Computation Cost : O(1)\n
 * It allocates the memory for queue and return allocated TaskQueue
 * @printEachElementFunc    : Call back function called for each data when print is called
 * @param heap              : Pointer to static heap
 *                          : NULL for dynamic heap
 * @return                  : Allocated TaskQueue (!!! Must be free using free) (OR) NULL if heap is full
 */
static TaskQueue *new(BuddyHeap *heap,void (*printEachElementFunc)(Task task)) {
    //Allocate memory for hash map
    TaskQueue *queue = allocateMemory(heap,sizeof(TaskQueue));

    //Invalid memory
    if (!validMemory(__func__, heap, queue))
        return NULL;

    //Heap is full
    if (queue == NULL)
        return NULL;
    queue->heap = heap;

    queue->printEachElement = printEachElementFunc;

    //Making both front and back null
    queue->front = NULL;
    queue->back = NULL;

    //Make initial size zero
    queue->size = 0;
    return queue;
}

/**
 * Computation Cost : O(1)\n
 * It adds the element to the end of the queue
 * @param queue     : TaskQueue
 * @param task      : Task to be added in queue
 * @return          : Same que (OR) NULL if heap is full or queue is null
 */
static TaskQueue *enqueue(TaskQueue *queue, Task task) {
    //If map is NULL then return NULL
    if (queue == NULL)
        return NULL;

    //Allocate Memory for newData
    TaskQueueData *newData = allocateMemory(queue->heap,sizeof(TaskQueueData));

    //If newData is invalid pointer return NULL
    if (!validMemory(__func__, queue->heap, newData))
        return NULL;

    //If heap is full then return NULL
    if (newData == NULL)
        return NULL;

    //Fill the task in data
    newData->task = task;
    newData->next = NULL;//make next data is empty

    //Get the last data
    if (queue->size == 0) {
        //If que is empty
        queue->front = newData;
        queue->back = newData;
    } else {
        queue->back->next = newData;
        queue->back = newData;
    }
    //If data is added then increase the size
    queue->size++;

    //Return same @queue
    return queue;
}

/**
 * Computation Cost : O(1)\n
 * It remove the element from the front of the queue and return it
 * @param queue     : TaskQueue
 * @return          : Element in front (OR) QUE_NULL if queue is empty or queue is null
 */
static Task dequeue(TaskQueue *queue) {
    //If map is NULL then return NULL
    if (queue == NULL)
        return TASK_NULL;

    //Get the last data
    if (queue->size == 0) {
        //If que is empty
        return TASK_NULL;
    } else {
        //Get the front data of queue
        TaskQueueData *frontData = queue->front;
        Task task = frontData->task;

        //Put front second data in the front
        queue->front = frontData->next;

        //Deallocate the allocated memory by front data
        freeMemory(queue->heap,frontData, sizeof(TaskQueueData));

        //Decrease the size of queue
        queue->size--;
        if (queue->size == 0) {
            queue->front = NULL;
            queue->back = NULL;
        }
        return task;
    }
}

/**
 * Computation Cost : O(1)\n
 * It returns the element from the front of the queue without removing it
 * @param queue     : TaskQueue
 * @return          : Element in front (OR) QUE_NULL if queue is empty or queue is null
 */
static Task peek(TaskQueue *queue) {
    //If map is NULL then return NULL
    if (queue == NULL)
        return TASK_NULL;

    if (queue->size == 0) {
        //If que is empty
        return TASK_NULL;
    } else {
        //Return the front element of queue
        return queue->front->task;
    }
}

/**
 * Computation Cost : O(n)\n
 * It delete all the data and free memories allocated by queue
 * @param queuePtr  : Address of pointer to queue
 * @return          : 1 for success (OR) 0 for failed
 */
static int freeQue(TaskQueue **queuePtr) {
    TaskQueue *queue = *queuePtr;
    //If queue is NULL
    if (queue == NULL)
        return 0;

    int size = queue->size;
    //If queue is empty
    if (size == 0) {
        //Free hash map memory
        freeMemory(queue->heap,queue, sizeof(TaskQueue));
        *queuePtr = NULL;
        return 1;
    }

    //Delete all data
    for (int i = 0; i < size; ++i)
        dequeue(queue);

    //Free memory for queue
    freeMemory(queue->heap,queue, sizeof(TaskQueue));

    *queuePtr = NULL;

    return 1; //freeing memory success
}

/**
 * This will print the contents of que
 * @param queue : TaskQueue to be printed
 */
static void print(TaskQueue *queue) {
    if (queue == NULL)
        return;

    if (queue->printEachElement == NULL)
        return;

    TaskQueueData *data = queue->front;
    for (int i = 0; i < queue->size; ++i) {
        if (data == NULL)
            break;
        queue->printEachElement(data->task);
        data = data->next;
    }
}

/**
 * This should be called regularly. It handles all the tasks and run it.
 * @param queue     : TaskQue
 */
static void execute(TaskQueue *queue){
    //If queue is NULL
    if (queue == NULL)
        return;

    int size = queue->size;
    //If queue is empty
    if (size == 0)
        return;

    Task task = dequeue(queue);
    if(task == TASK_NULL)//If task is corrupt
        return;

    //If task is not corrupt run it
    TaskStatus status = task();

    //If task is busy then again push the task at the end
    if(status==TASK_BUSY)
        enqueue(queue, task);
}

/**
 * This return allocated memory for queue till now
 * @return  : Allocated memories
 */
static int getAllocatedMemories() {
    return allocatedMemory;
}
struct TaskQueControl StaticTaskQueue = {
        .new=new,
        .enqueue = enqueue,
        .dequeue = dequeue,
        .peek = peek,
        .free=freeQue,
        .print=print,
        .execute = execute,
        .getAllocatedMemories=getAllocatedMemories
};



