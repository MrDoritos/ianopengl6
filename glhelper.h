#pragma once 
#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <GLFW/glfw3.h>
#include <GL/glew.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "vertex.h"

inline std::string* readShaderFile(const char* filePath) {
    std::ifstream shaderFile(filePath);
    std::stringstream shaderStream;
    shaderStream << shaderFile.rdbuf();
    std::string *ret = new std::string(shaderStream.str());
    shaderFile.close();
    return ret;
}

bool compileShader(GLuint shaderId, const char* shaderSource) {
    glShaderSource(shaderId, 1, &shaderSource, nullptr);
    glCompileShader(shaderId);

    GLint compileStatus;
    glGetShaderiv(shaderId, GL_COMPILE_STATUS, &compileStatus);
    if (compileStatus == GL_FALSE) {
        GLint logLength;
        glGetShaderiv(shaderId, GL_INFO_LOG_LENGTH, &logLength);
        std::string errorMessage(logLength, ' ');
        glGetShaderInfoLog(shaderId, logLength, nullptr, &errorMessage[0]);
        std::cout << "Shader compilation error: " << errorMessage << std::endl;
        return false;
    }

    return true;
}

bool linkProgram(GLuint programId, GLuint vertexShaderId, GLuint fragmentShaderId) {
    glAttachShader(programId, vertexShaderId);
    glAttachShader(programId, fragmentShaderId);
    glLinkProgram(programId);

    GLint linkStatus;
    glGetProgramiv(programId, GL_LINK_STATUS, &linkStatus);
    if (linkStatus == GL_FALSE) {
        GLint logLength;
        glGetProgramiv(programId, GL_INFO_LOG_LENGTH, &logLength);
        std::string errorMessage(logLength, ' ');
        glGetProgramInfoLog(programId, logLength, nullptr, &errorMessage[0]);
        std::cout << "Program linking error: " << errorMessage << std::endl;
        return false;
    }

    glDetachShader(programId, vertexShaderId);
    glDetachShader(programId, fragmentShaderId);

    return true;
}

GLuint figureOutShaders(const char* prefix) {
        // Load and compile vertex shader
    GLuint program_ret;

    std::string vertex_shader_path = std::string(prefix) + std::string("vertex_shader.glsl");
    std::string fragment_shader_path = std::string(prefix) + std::string("fragment_shader.glsl");

    std::string *vertexShaderSource =(readShaderFile(vertex_shader_path.c_str()));
    const char* vertexShaderSourcePtr = vertexShaderSource->c_str();
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    if (!compileShader(vertexShader, vertexShaderSourcePtr)) {
        glfwTerminate();
        return 0;
    }

    // Load and compile fragment shader
    std::string *fragmentShaderSource = (readShaderFile(fragment_shader_path.c_str()));
    const char* fragmentShaderSourcePtr = fragmentShaderSource->c_str();
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    if (!compileShader(fragmentShader, fragmentShaderSourcePtr)) {
        glfwTerminate();
        return 0;
    }

    // Create shader program and link shaders
    program_ret = glCreateProgram();
    if (!linkProgram(program_ret, vertexShader, fragmentShader)) {
        glfwTerminate();
        return 0;
    }

    // Use the shader program
    glUseProgram(program_ret);

    return program_ret;
}

int loadTexture(const char* path) {
    // Load texture using stb_image
    int width, height, channels;
    unsigned char* image = stbi_load(path, &width, &height, &channels, 0);
    if (!image) {
        std::cout << "Failed to load texture" << std::endl;
        glfwTerminate();
        return -1;
    }

    // Print texture information
    std::cout << "Loaded texture: width = " << width << ", height = " << height << ", channels = " << channels << std::endl;


        // Load texture into OpenGL
    GLuint textureId;
    glGenTextures(1, &textureId);
    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
    glGenerateMipmap(GL_TEXTURE_2D);

    // Set texture wrapping and filtering options
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    return textureId;
}
