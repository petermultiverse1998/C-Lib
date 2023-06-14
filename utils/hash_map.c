//
// Created by peter on 6/13/2023.
//

#include "hash_map.h"
#include <stdint.h>
#include "stdio.h"
#include "stdlib.h"

#define MAP_SIZE HASH_MAP_SIZE
#define MAX_LOOP 1000

/**
 * This converts key into hasKey
 * @param key   : Actual Key
 * @return      : Hash Key
 */
static int hashFunc(int key) {
    return key % MAP_SIZE;
}

/**
 * It allocates the memory and return pointer to it
 * @param sideInByte    : Size in bytes
 * @return              : Pointer to allocated memory
 *                      : NULL if there exist no memory for allocation
 */
static void *allocateMemory(int sideInByte) {
//    return NULL;
    return malloc(sideInByte);
}

/**
 * It free the allocated memory
 * @param pointer       : Pointer to allocated Memory
 * @param sizeInByte    : Size to be freed
 * @return              : 1 for success (OR) 0 for failed
 */
static int freeMemory(void *pointer, int sizeInByte) {
    (void)sizeInByte;
//    if(sizeInByte==0)
//        return 1;
    free(pointer);
    return 1;
}

struct Entry {
    int key;
    Type value;
    struct Entry *nextEntryPtr;
};
typedef struct Entry Entry;

struct HashMap {
    int size;
    Entry *entries[MAP_SIZE];
};

/**
 * Computation Cost : O(1)\n
 * It allocates the memory for HashMap and return allocated HashMap
 * @return : Allocated HashMap (!!! Must be free using free) (OR) NULL if heap is full
 */
static HashMap *new() {
    //Allocate memory for hash map
    HashMap *map = allocateMemory(sizeof(HashMap));

    //Heap is full
    if (map == NULL)
        return NULL;

    //Filling all entries with NULL
    for (int i = 0; i < MAP_SIZE; ++i)
        map->entries[i] = NULL;

    //Make initial size zero
    map->size = 0;
    return map;
}

/**
 * Prints the contents of hash map
 * @param map  : Hash map
 */
static void print(HashMap *map) {
    if (map == NULL) {
        printf("map is NULL!!!\n");
        return;
    }

    for (int i = 0; i < MAP_SIZE; ++i) {
        printf("%4d: ", i);
        Entry *entry = map->entries[i];
        for (int j = 0; j < MAX_LOOP; ++j) {
            if (entry == NULL)
                break;
            printf("0x%02x >> ", entry->key);
            entry = entry->nextEntryPtr;
        }
        printf("\n");
    }
}

/**
 * Computation Cost : O(1)\n
 * It inserts the value in map
 * @param map   : HashMap
 * @param key   : Key for value
 * @param value : Value to be inserted of type Type
 * @return      : Same map (OR) NULL if heap is full or map is null
 */
static HashMap *insert(HashMap *map, int key, Type value) {
    //If map is NULL then return NULL
    if (map == NULL)
        return NULL;

    //Allocate Memory for @Entry
    Entry *newEntry = allocateMemory(sizeof(Entry));

    //If heap is full then return NULL
    if (newEntry == NULL)
        return NULL;

    //Fill the data in entry
    newEntry->key = key;
    newEntry->value = value;
    newEntry->nextEntryPtr = NULL;//make next entry is empty

    //Calculate hashKey
    int hashKey = hashFunc(key);

    //Get the top of entry of corresponding hashmap
    Entry *entry = map->entries[hashKey];

    if (entry == NULL) {
        //If top entry is empty fill the entry
        map->entries[hashKey] = newEntry;
    } else {
        //Check if key already exist
        if (entry->key == key) {
            //Only update value
            entry->value = value;
            return map;
        }

        //If top entry is not empty get next entry
        Entry *nextEntry = entry->nextEntryPtr;
        int l = 0;
        for (; l < MAX_LOOP; ++l) {
            //If this entry is empty(@NULL) break
            if (nextEntry == NULL) {
                entry->nextEntryPtr = newEntry;
                break;
            }
            //If this entry is not empty select next entry
            entry = nextEntry;
            nextEntry = nextEntry->nextEntryPtr;
        }

        //If @MAX_LOOP exceeds then entry is not inserted
        if (l >= MAX_LOOP) {
            //Free allocated memory for @newEntry
            free(newEntry);
            return NULL;
        }

        //If entry is inserted then increase the size
        map->size++;
    }
    //Return same @map
    return map;
}

/**
 * Computation Cost : O(1)\n
 * It searches the value corresponding to the key and return value
 * @param map   : HashMap
 * @param key   : Key for value
 * @return      : Value corresponding to the key (OR) NULL if key doesn't exist found
 */
static Type get(HashMap *map, int key) {
    //If map is NULL return NULL
    if (map == NULL)
        return map;

    //Calculate hashKey
    int hashKey = hashFunc(key);

    //Get the top of entry for corresponding hashmap
    Entry *entry = map->entries[hashKey];
    for (int l = 0; l < MAX_LOOP; ++l) {
        //If this entry is empty value doesn't exist in map so return NULL
        if (entry == NULL)
            return NULL;

        //If key is found return value
        if (entry->key == key)
            return entry->value;

        //If key is not found then go to next entry
        entry = entry->nextEntryPtr;
    }
    //If loop exceeds the limit and no key is found then return NULL
    return NULL;
}

/**
 * Computation Cost : O(1)\n
 * It deletes the values in hash map corresponding to given key
 * @param map   : HashMap
 * @param key   : Key for value
 * @return      : Same map (OR) NULL if key is not found or map is NULL
 */
static HashMap *delete(HashMap *map, int key) {
    //If map is NULL return NULL
    if (map == NULL)
        return NULL;

    //Calculate hashKey
    int hashKey = hashFunc(key);

    //Get the top of entry of corresponding hashmap
    Entry *entry = map->entries[hashKey];
    //Initially previous entry is NULL
    Entry *prevEntry = NULL;
    for (int l = 0; l < MAX_LOOP; ++l) {
        //If entry is empty key doesn't exist
        if (entry == NULL)
            break;

        //Key found
        if (entry->key == key) {
            //Copy next entry
            Entry *nextEntry = entry->nextEntryPtr;

            if (prevEntry == NULL)
                // If prevEntry is NULL then this entry is top entry
                map->entries[hashKey] = nextEntry;
            else
                // If prevEntry is not NULL then this entry is branch entry
                prevEntry->nextEntryPtr = nextEntry;

            //Decrease the size
            map->size--;
            return map;
        }

        //Go to next entry if key is not found
        prevEntry = entry;
        entry = entry->nextEntryPtr;
    }
    //If key is not found return NULL
    return NULL;
}

/**
 * Computation Cost : O(n)\n
 * It gives all the keys in the map
 * @param map   : HashMap
 * @param keys  : Array of integer of size equal to size of map
 * @return      : same map (OR) NULL if map is NULL
 */
static HashMap *getKeys(HashMap *map, int keys[]) {
    //If map is NULL return NULL
    if (map == NULL)
        return NULL;

    int size = map->size;
    int keyIndex = 0;
    for (int entryIndex = 0; entryIndex < size; ++entryIndex) {
        //Get top entry
        Entry *entry = map->entries[entryIndex];

        for (int l = 0; l < MAX_LOOP; ++l) {
            //If entry is empty then break loop
            if (entry == NULL)
                break;

            //If entry is not empty
            keys[keyIndex++] = entry->key;// Copy key and increase keyIndex
            entry = entry->nextEntryPtr;//Go to next entry
        }
    }
    return map;
}

/**
 * Computation Cost : O(1)\n
 * It searches if key exist or not in map
 * @param map   : HashMap
 * @param key   : Key to be searched
 * @return      : 1 if key exist (OR) 0 if key doesn't exist
 */
static int isKeyExist(HashMap *map, int key) {
    //If map is NULL return false
    if (map == NULL)
        return 0;

    //Calculate hashKey
    int hashKey = hashFunc(key);

    //Get the top of entry for corresponding hashmap
    Entry *entry = map->entries[hashKey];
    for (int l = 0; l < MAX_LOOP; ++l) {
        //If this entry is empty key doesn't exist in map
        if (entry == NULL)
            return 0;

        //If key is found return true
        if (entry->key == key)
            return 1;

        //If key is not found then go to next entry
        entry = entry->nextEntryPtr;
    }
    //If loop exceeds the limit and no key is found then return false
    return 0;
}

/**
 * Computation Cost : O(n^2)\n
 * It delete all the entries and free memories allocated by map
 * @param map   : HashMap
 * @return      : 1 for success (OR) 0 for failed
 */
static int freeMap(HashMap *map) {
    //If map is NULL
    if (map == NULL)
        return 0;

    int size = map->size;
    //If map is empty
    if (size == 0) {
        //Free hash map memory
        freeMemory(map, sizeof(HashMap));
        return 1;
    }

    //If map is not empty
    int keys[size];

    //Get all keys
    getKeys(map, keys);


    //Delete all entries
    for (int i = 0; i < size; ++i)
        delete(map, keys[i]);

    //Free memory for hash map
    freeMemory(map, sizeof(HashMap));

    return 1; //freeing memory success
}


struct HashMapControl StaticHashMap = {
        .new=new,
        .insert=insert,
        .get=get,
        .delete=delete,
        .getKeys=getKeys,
        .isKeyExist=isKeyExist,
        .free=freeMap,
        .print=print
};