#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "include/stb_image.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>

#include "include/Shader.h"
#include "include/Camera.h"
#include "include/Mesh.h"
//#include "include/Model.h"

#include <iostream>
#include <string>
#include <cmath>

#include "Simulation.h"
#include "particle.h"

// settings
const unsigned int SCR_WIDTH = 1200;
const unsigned int SCR_HEIGHT = 1200;

// camera
Camera camera(glm::vec3(5.0f, 2.0f, 9.0f));
float lastX = (float)SCR_WIDTH / 2.0;
float lastY = (float)SCR_HEIGHT / 2.0;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);

const int FPS = 30;
std::string str_cmd = "ffmpeg -r " + std::to_string(FPS) + " -f rawvideo -pix_fmt rgba -s "
+ std::to_string(SCR_WIDTH) + "x" + std::to_string(SCR_HEIGHT)
+ " -i - -threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip ./3D.mp4";
const char* cmd = str_cmd.c_str();
FILE* ffmpeg = _popen(cmd, "wb");
int* buffer = new int[SCR_WIDTH * SCR_HEIGHT];

//GLFWwindow* glfwSetup()
//{
//    // glfw: initialize and configure
//// ------------------------------
//    glfwInit();
//    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
//    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
//    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
//
//    // glfw window creation
//    // --------------------
//    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
//    if (window == NULL)
//    {
//        std::cout << "Failed to create GLFW window" << std::endl;
//        glfwTerminate();
//        //return -1;
//    }
//    glfwMakeContextCurrent(window);
//    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
//    glfwSetCursorPosCallback(window, mouse_callback);
//    glfwSetScrollCallback(window, scroll_callback);
//
//    // tell GLFW to capture our mouse
//    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
//
//    // tell GLFW to capture our mouse
//    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
//
//    return window;
//}

std::shared_ptr<Mesh> get_sphere(float r, size_t slices, size_t stacks) {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    const float kPi = std::atan(1.0f) * 4;

    float phi_step = kPi * 2 / slices;
    float theta_step = kPi / stacks;

    for (size_t vi = 0; vi < stacks; vi++) {  // vertical loop
        float theta = vi * theta_step;
        float z = r * cosf(theta);
        for (size_t hi = 0; hi < slices; hi++) {  // horizontal loop
            float phi = hi * phi_step;
            glm::vec3 p(r * cosf(phi) * sinf(theta), r * sinf(phi) * sinf(theta), z);
            Vertex vertex;
            vertex.Position = p;
            vertex.Normal = glm::normalize(p);
            vertices.push_back(vertex);
        }
    }

    for (size_t vi = 0; vi < stacks; vi++)
        for (size_t hi = 0; hi < slices; hi++) {
            auto t1 = (unsigned int)(vi * slices + hi);
            auto t2 = (unsigned int)(vi * slices + hi + 1);
            auto t3 = (unsigned int)((vi + 1) * slices + hi + 1);
            auto t4 = (unsigned int)((vi + 1) * slices + hi);
            indices.insert(indices.end(), { t1, t2, t3 });
            indices.insert(indices.end(), { t1, t3, t4 });
        }

    std::shared_ptr<Mesh> sphere_mesh = std::make_unique<Mesh>(vertices, indices);

    return sphere_mesh;
}

//std::unique_ptr<Mesh> get_square() {
//    std::vector<Vertex> vertices;
//    std::vector<unsigned int> indices;
//
//    Vertex vertex1;
//    Vertex vertex2;
//    Vertex vertex3;
//    Vertex vertex4;
//    vertex1.Position = glm::vec3(-0.05f, 0.05f, 0.0f);
//    vertex2.Position = glm::vec3(0.05f, -0.05f, 0.0f);
//    vertex3.Position = glm::vec3(-0.05f, -0.05f, 0.0f);
//    vertex4.Position = glm::vec3(0.05f, 0.05f, 0.0f);
//    vertex1.Normal = glm::vec3(0.0f, 0.0f, 1.0f);
//    vertex2.Normal = glm::vec3(0.0f, 0.0f, 1.0f);
//    vertex3.Normal = glm::vec3(0.0f, 0.0f, 1.0f);
//    vertex4.Normal = glm::vec3(0.0f, 0.0f, 1.0f);
//    vertices.push_back(vertex1);
//    vertices.push_back(vertex2);
//    vertices.push_back(vertex3);
//    vertices.push_back(vertex4);
//    indices.insert(indices.end(), { 0, 1, 2 });
//    indices.insert(indices.end(), { 0, 1, 3 });
//
//    std::unique_ptr<Mesh> square_mesh = std::make_unique<Mesh>(vertices, indices);
//
//    return square_mesh;
//}

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile shaders
    // -------------------------
    Shader shader("./shader/shader.vs", "./shader/shader.fs");
    Shader bordershader("./shader/shader_ez.vs", "./shader/shader_ez.fs");

    // generate a list of 100 quad locations/translation-vectors
    // ---------------------------------------------------------

    //glm::mat4* modelMatrices;
    //modelMatrices = new glm::mat4[100];
    //int index = 0;
    //float offset = 0.1f;
    //for (int y = -10; y < 10; y += 2)
    //{
    //    for (int x = -10; x < 10; x += 2)
    //    {
    //        glm::vec3 translation;
    //        translation.x = (float)x / 10.0f + offset;
    //        translation.y = (float)y / 10.0f + offset;
    //        translation.z = 0.0f; //cannot be greater than 0.1f

    //        glm::mat4 model = glm::mat4(1.0f);
    //        //std::cout << glm::to_string(model) << std::endl;

    //        model = glm::translate(model, translation);
    //        //std::cout << glm::to_string(model) << std::endl;

    //        modelMatrices[index++] = model;
    //    }
    //}

    //auto sphere = get_square();

    //unsigned int instanceVBO;
    //glGenBuffers(1, &instanceVBO);
    //glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(glm::mat4) * 100, &modelMatrices[0], GL_STATIC_DRAW);

    //glBindVertexArray(sphere->VAO);

    //glEnableVertexAttribArray(3);
    //glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)0);
    //glEnableVertexAttribArray(4);
    //glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)(sizeof(glm::vec4)));
    //glEnableVertexAttribArray(5);
    //glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)(2 * sizeof(glm::vec4)));
    //glEnableVertexAttribArray(6);
    //glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)(3 * sizeof(glm::vec4)));

    //glVertexAttribDivisor(3, 1);
    //glVertexAttribDivisor(4, 1);
    //glVertexAttribDivisor(5, 1);
    //glVertexAttribDivisor(6, 1);

    //glBindVertexArray(0);


    //glm::vec3 a(1.0f, 2.0f, 3.0f);
    //glm::vec3 b(2.0f, 2.0f, 2.0f);
    //glm::vec3 c = a * b;
    //float d = glm::dot(a,b);
    //glm::mat3 e = glm::outerProduct(a, b);
    //std::cout << glm::to_string(c) << std::endl;
    //std::cout << d << std::endl;
    //std::cout << glm::to_string(e) << std::endl;

    //==========================================================================================
    auto sphere = get_sphere(0.02f, 20, 20);

    glm::vec3 v0(0.0f,0.0f,0.0f);
    glm::mat3 B0(0.0f);
    float mass = 0.0005f;
    float V0 = 1.14f;

    float dt = 0.001f;

    std::vector<Water> outParticles;

    srand((unsigned)time(NULL));
    for (int p = 0; p < 10000; p++)
    {
        float r1 = ((float)rand() / (RAND_MAX));			// random number
        float r2 = ((float)rand() / (RAND_MAX));			// random number
        float r3 = ((float)rand() / (RAND_MAX));
        float theta = 2.0f * 3.14159265f * r1;
        float phi = 3.14159265f * r2;
        float radius = 1.0f * r3;
        float x = radius * sin(phi) * cos(theta);
        float z = radius * sin(phi) * sin(theta);
        float y = radius * cos(phi);
        glm::vec3 pos(50.0f + x, 5.0f + y, 50.0f + z);
        Water the_water(pos, v0, B0, mass, V0);
        outParticles.push_back(the_water);
    }

    Simulation simulation(outParticles, dt);

    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    glm::mat4 view = camera.GetViewMatrix();

    //==========================================================================================


    // render loop
    // -----------
    int count = 0;
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
// --------------------
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();

        // draw 100 instanced quads
        //shader.use();
        //shader.setMat4("projection", projection);
        //shader.setMat4("view", view);
        ////glBindVertexArray(quadVAO);
        //glBindVertexArray(sphere->VAO);
        ////glDrawArraysInstanced(GL_TRIANGLES, 0, 6, 100); // 100 triangles of 6 vertices each
        //glDrawElementsInstanced(GL_TRIANGLES, sphere->indices.size(), GL_UNSIGNED_INT, 0, 100);
        //glBindVertexArray(0);
        //srand((unsigned)time(NULL));
        //if (count <= 21000 && count >=20900) {
        //    if (count % 5 == 0) {
        //        for (int i = 0; i < 50; i++) {
        //            double r1 = ((double)rand() / (RAND_MAX));
        //            double r2 = ((double)rand() / (RAND_MAX));
        //            float theta = 2.0f * 3.14159265f * r1;
        //            float radius = 5.0f * r2;
        //            float x = 50.0f + radius * cos(theta);
        //            float z = 50.0f + radius * sin(theta);
        //            glm::vec3 pos(x, 5.0f, z);
        //            Water the_water(pos, v0, B0, mass, V0);
        //            simulation.add_particle(the_water);
        //        }
        //    }
        //}

        simulation.forward_one_step();

        if (count % 30 == 0) {
            simulation.render_particles(shader, sphere, projection, view);
            simulation.render_borders(bordershader, projection, view);
            glfwSwapBuffers(window);

            glReadPixels(0, 0, SCR_WIDTH, SCR_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
            fwrite(buffer, sizeof(int) * SCR_WIDTH * SCR_HEIGHT, 1, ffmpeg);

            glfwPollEvents();
        }
        simulation.reset_grid();

        count++;

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        //glfwSwapBuffers(window);
        //glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    //glDeleteVertexArrays(1, &quadVAO);
    //glDeleteBuffers(1, &quadVBO);

    _pclose(ffmpeg);
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}