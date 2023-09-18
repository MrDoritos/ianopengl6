#pragma once
#include <vector>

#include "vertex.h"

struct block_t;
struct fullblock_t;
struct blockstate_t;

struct blockstate_t {
    blockstate_t(short blockId) {
        this->blockId = blockId;
    }
    blockstate_t() {
    }
    short blockId;
};

struct block_t {
    static std::vector<block_t*> blocks;

    static block_t *air;
    static block_t *stone;
    static block_t *grass;
    static block_t *dirt;

    block_t() {
        blockId = 0;
        
    }

    block_t(int blockId, char t0, char t1, char t2, char t3) {
        this->blockId = blockId;
        textureAtlas[0] = t0;
        textureAtlas[1] = t1;
        textureAtlas[2] = t2;
        textureAtlas[3] = t3;
        blocks.push_back(this);
    }


    virtual blockstate_t getDefaultState() {
        blockstate_t ret;
        ret.blockId = blockId;
        return ret;
    }

    virtual unsigned char getFaceBrightness(int face, fullblock_t fb);

	virtual unsigned char getTexture(int triangle, int vertex, int index) {						
		return textureAtlas[map[triangle % 2][vertex][index]];
	}

    int blockId;  
	unsigned char textureAtlas[4];
};

std::vector<block_t*> block_t::blocks;
block_t *block_t::air;
block_t *block_t::grass;
block_t *block_t::stone;
block_t *block_t::dirt;

blockstate_t nothing(-1);

struct fullblock_t {
    fullblock_t() {
        state = &nothing;
        block = block_t::stone;
    }
    vec3d position;
    blockstate_t *state;
    block_t *block;
};

unsigned char block_t::getFaceBrightness(int face, fullblock_t fb) {
    switch (face) {
        case TOP_FACE:
            return 255;
        case BOTTOM_FACE:
            return 127;
        default:
            return 191;
    }
}
