#include <imgui.h>
#include <implot.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif

struct TriangleInspectData
{
    double plot_x[3];
    double plot_y[3];
    double lengths[3];
    double angles[3];
    double area;

    double side_equations[3][3];

    double heights[3];
    double heights_x[3];
    double heights_y[3];
    double height_equations[3][3];

    double medians_x[3];
    double medians_y[3];
    double median_lengths[3];

    double centroid_x;
    double centroid_y;
};

static double Vector2Length(double x, double y)
{
    return sqrt((x * x) + (y * y));
}

static double Vector2DotProduct(double x0, double y0, double x1, double y1)
{
    return (x0 * x1) + (y0 * y1);
}

static void VectorExtend(double bx, double by, double x, double y, double scale_x, double scale_y, double* ox, double* oy)
{
    double len = Vector2Length(y, x);
    double dir_norm[2] = { x / len, y / len };
    *ox = bx + (dir_norm[0] * scale_x);
    *oy = by + (dir_norm[1] * scale_y);
}

static void BuildAndSimplifyEquation(double a, double b, double c, double* eq)
{
    eq[0] = a;
    eq[1] = b;
    eq[2] = c;

    double r = remainder(eq[0], eq[1]);
    if (r == 0.0)
    {
        const double d = eq[1];
        for (int j = 0; j < 3; j++)
            eq[j] /= d;
    }
    else
    {
        r = remainder(eq[1], eq[0]);
        if (r == 0.0)
        {
            const double d = eq[0];
            for (int j = 0; j < 3; j++)
                eq[j] /= d;
        }
    }

    if (eq[0] < 0.0)
    {
        for (int j = 0; j < 3; j++)
            eq[j] = -eq[j];
    }
}

static void FormatEquation(char* buf, size_t buf_count, double* eq)
{
    const double abs_a = fabs(eq[0]);
    const double abs_b = fabs(eq[1]);

    char fmt[32];
    if (eq[2] == 0.0)
    {
        const char* y_fmt;
        if (abs_b == 1.0)
            y_fmt = eq[1] < 0.0 ? "- " : "+ ";
        else if (abs_b != 0.0 && eq[1] < 0.0)
            y_fmt = "- %g";
        else
            y_fmt = "+ %g";

        snprintf(fmt, sizeof(fmt), "%sx %sy = 0",
            abs_a != 1.0 ? "%g" : "",
            y_fmt
        );

        if (abs_a != 1.0)
        {
            if (abs_b != 1.0)
                snprintf(buf, buf_count, fmt, fabs(eq[0]), fabs(eq[1]));
            else
                snprintf(buf, buf_count, fmt, fabs(eq[0]));
        }
        else
        {
            if (abs_b != 1.0)
                snprintf(buf, buf_count, fmt, fabs(eq[1]));
            else
                snprintf(buf, buf_count, fmt);
        }
    }
    else
    {
        const char* y_fmt;
        if (abs_b == 1.0)
            y_fmt = eq[1] < 0.0 ? "- " : "+ ";
        else if (abs_b != 0.0 && eq[1] < 0.0)
            y_fmt = "- %g";
        else
            y_fmt = "+ %g";

        snprintf(fmt, sizeof(fmt), "%sx %sy %s = 0",
            abs_a != 1.0 ? "%g" : "",
            y_fmt,
            eq[2] < 0.0 ? "- %g" : "+ %g"
        );

        if (abs_a != 1.0)
        {
            if (abs_b != 1.0)
                snprintf(buf, buf_count, fmt, fabs(eq[0]), fabs(eq[1]), fabs(eq[2]));
            else
                snprintf(buf, buf_count, fmt, fabs(eq[0]), fabs(eq[2]));
        }
        else
        {
            if (eq[1] != 1.0)
                snprintf(buf, buf_count, fmt, fabs(eq[1]), fabs(eq[2]));
            else
                snprintf(buf, buf_count, fmt, fabs(eq[2]));
        }
    }
}

static void UpdateTriangleInspectData(double x0, double y0, double x1, double y1, double x2, double y2, struct TriangleInspectData* data)
{
    data->plot_x[0] = x0;
    data->plot_x[1] = x1;
    data->plot_x[2] = x2;

    data->plot_y[0] = y0;
    data->plot_y[1] = y1;
    data->plot_y[2] = y2;

    data->lengths[0] = Vector2Length(data->plot_x[1] - data->plot_x[0], data->plot_y[1] - data->plot_y[0]); // AB
    data->lengths[1] = Vector2Length(data->plot_x[2] - data->plot_x[1], data->plot_y[2] - data->plot_y[1]); // BC
    data->lengths[2] = Vector2Length(data->plot_x[0] - data->plot_x[2], data->plot_y[0] - data->plot_y[2]); // AC

    // law of cosines
    data->angles[0] = acos(((data->lengths[0] * data->lengths[0]) + (data->lengths[2] * data->lengths[2]) - (data->lengths[1] * data->lengths[1]))
        / (2.0 * data->lengths[0] * data->lengths[2])); // A
    data->angles[1] = acos(((data->lengths[0] * data->lengths[0]) + (data->lengths[1] * data->lengths[1]) - (data->lengths[2] * data->lengths[2]))
        / (2.0 * data->lengths[0] * data->lengths[1])); // B
    data->angles[2] = acos(((data->lengths[2] * data->lengths[2]) + (data->lengths[1] * data->lengths[1]) - (data->lengths[0] * data->lengths[0]))
        / (2.0 * data->lengths[2] * data->lengths[1])); // C

    double area = ((data->plot_x[1] - data->plot_x[0]) * (data->plot_y[2] - data->plot_y[0]) - (data->plot_x[2] - data->plot_x[0]) * (data->plot_y[1] - data->plot_y[0]));
    area = (area < 0 ? -area : area) / 2.0;
    data->area = area;

    double tmp0, tmp1;
    tmp0 = -(x1 - x0);
    tmp1 = y1 - y0;

    BuildAndSimplifyEquation(tmp1, tmp0, -((tmp0 * y0) + (tmp1 * x0)), data->side_equations[0]);

    tmp0 = -(x2 - x0);
    tmp1 = y2 - y0;

    BuildAndSimplifyEquation(tmp1, tmp0, -((tmp0 * y0) + (tmp1 * x0)), data->side_equations[1]);

    tmp0 = -(x2 - x1);
    tmp1 = y2 - y1;

    BuildAndSimplifyEquation(tmp1, tmp0, -((tmp0 * y1) + (tmp1 * x1)), data->side_equations[2]);

    data->heights[0] = (2.0 * data->area) / data->lengths[1]; // AH1
    data->heights[1] = (2.0 * data->area) / data->lengths[2]; // BH2
    data->heights[2] = (2.0 * data->area) / data->lengths[0]; // CH2

    tmp0 = -(y2 - y1);
    tmp1 = x1 - x2;

    BuildAndSimplifyEquation(tmp1, tmp0, -((tmp0 * y0) + (tmp1 * x0)), data->height_equations[0]);

    tmp0 = -(y2 - y0);
    tmp1 = x0 - x2;

    BuildAndSimplifyEquation(tmp1, tmp0, -((tmp0 * y1) + (tmp1 * x1)), data->height_equations[1]);

    tmp0 = -(y1 - y0);
    tmp1 = x0 - x1;

    BuildAndSimplifyEquation(tmp1, tmp0, -((tmp0 * y2) + (tmp1 * x2)), data->height_equations[2]);

    // height bases

    // proj BC(BA)
    tmp0 = Vector2DotProduct(x0 - x1, y0 - y1, x2 - x1, y2 - y1); // BA * BC
    tmp1 = Vector2DotProduct(x2 - x1, y2 - y1, x2 - x1, y2 - y1); // BC * BC

    data->heights_x[0] = ((x0 * tmp1) - (((x0 - x1) * tmp1) - (tmp0 * (x2 - x1)))) / tmp1;
    data->heights_y[0] = ((y0 * tmp1) - (((y0 - y1) * tmp1) - (tmp0 * (y2 - y1)))) / tmp1;

    // proj AC(AB)
    tmp0 = Vector2DotProduct(x1 - x0, y1 - y0, x2 - x0, y2 - y0); // AB * AC
    tmp1 = Vector2DotProduct(x2 - x0, y2 - y0, x2 - x0, y2 - y0); // AC * AC

    data->heights_x[1] = ((x1 * tmp1) - (((x1 - x0) * tmp1) - (tmp0 * (x2 - x0)))) / tmp1;
    data->heights_y[1] = ((y1 * tmp1) - (((y1 - y0) * tmp1) - (tmp0 * (y2 - y0)))) / tmp1;

    // proj AB(AC)
    tmp0 = Vector2DotProduct(x2 - x0, y2 - y0, x1 - x0, y1 - y0); // AC * AB
    tmp1 = Vector2DotProduct(x1 - x0, y1 - y0, x1 - x0, y1 - y0); // AB * AB

    data->heights_x[2] = ((x2 * tmp1) - (((x2 - x0) * tmp1) - (tmp0 * (x1 - x0)))) / tmp1;
    data->heights_y[2] = ((y2 * tmp1) - (((y2 - y0) * tmp1) - (tmp0 * (y1 - y0)))) / tmp1;

    // M1
    data->medians_x[0] = (x1 + x2) / 2.0;
    data->medians_y[0] = (y1 + y2) / 2.0;

    // M2
    data->medians_x[1] = (x0 + x2) / 2.0;
    data->medians_y[1] = (y0 + y2) / 2.0;

    // M3
    data->medians_x[2] = (x0 + x1) / 2.0;
    data->medians_y[2] = (y0 + y1) / 2.0;

    tmp0 = data->medians_x[0] - x0;
    tmp1 = data->medians_y[0] - y0;
    data->median_lengths[0] = sqrt((tmp0 * tmp0) + (tmp1 * tmp1));

    tmp0 = data->medians_x[1] - x1;
    tmp1 = data->medians_y[1] - y1;
    data->median_lengths[1] = sqrt((tmp0 * tmp0) + (tmp1 * tmp1));

    tmp0 = data->medians_x[2] - x2;
    tmp1 = data->medians_y[2] - y2;
    data->median_lengths[2] = sqrt((tmp0 * tmp0) + (tmp1 * tmp1));

    data->centroid_x = (x0 + x1 + x2) / 3.0;
    data->centroid_y = (y0 + y1 + y2) / 3.0;
}

static void RenderTriangleWindow()
{
    static float radius = 25.f;
    ImGui::SliderFloat("Angle Distance", &radius, 0.f, 100.f, "%.3f%%");

    static bool filled_angles = false;
    ImGui::Checkbox("Filled Angles", &filled_angles);

    static bool show_medians = false;
    ImGui::Checkbox("Show Medians", &show_medians);

    static bool show_heights = false;
    ImGui::Checkbox("Show Heights", &show_heights);

    static bool show_lengths = false;
    ImGui::Checkbox("Show Side Lenghts", &show_lengths);

    static double mypoints[3][2] = {
        { 2.0, -4.0 },
        { -2.0, -2.0 },
        { 4.0, 1.0 }
    };

    static struct TriangleInspectData data;
    static bool wants_data_update = true;

    wants_data_update = ImGui::InputDouble("x0", &mypoints[0][0]) || wants_data_update;
    wants_data_update = ImGui::InputDouble("y0", &mypoints[0][1]) || wants_data_update;
    
    wants_data_update = ImGui::InputDouble("x1", &mypoints[1][0]) || wants_data_update;
    wants_data_update = ImGui::InputDouble("y1", &mypoints[1][1]) || wants_data_update;

    wants_data_update = ImGui::InputDouble("x2", &mypoints[2][0]) || wants_data_update;
    wants_data_update = ImGui::InputDouble("y2", &mypoints[2][1]) || wants_data_update;

    if (wants_data_update)
    {
        UpdateTriangleInspectData(mypoints[0][0], mypoints[0][1], mypoints[1][0], mypoints[1][1], mypoints[2][0], mypoints[2][1], &data);
        wants_data_update = false;
    }

    if (ImPlot::BeginPlot("My plot"))
    {
        ImPlot::PushPlotClipRect();
        ImDrawList* draw_list = ImPlot::GetPlotDrawList();

        for (int i = 0; i < 3; i++)
        {
            const int ni = i == 2 ? 0 : i + 1;
            const int pi = i == 0 ? 2 : i - 1;

            // angles
            const double bx = data.plot_x[i];
            const double by = data.plot_y[i];
            const ImVec2 bp = ImPlot::PlotToPixels({ bx, by });

            const double nx = data.plot_x[ni];
            const double ny = data.plot_y[ni];
            const ImVec2 np = ImPlot::PlotToPixels({ nx, ny });

            const double px = data.plot_x[pi];
            const double py = data.plot_y[pi];
            const ImVec2 pp = ImPlot::PlotToPixels({ px, py });

            const ImVec2 hp = ImPlot::PlotToPixels({ data.heights_x[i], data.heights_y[i] });

            double np_len = Vector2Length((np.x - bp.x) * 0.5f, (np.y - bp.y) * 0.5f);
            double pp_len = Vector2Length((pp.x - bp.x) * 0.5f, (pp.y - bp.y) * 0.5f);
            double hp_len = Vector2Length((hp.x - bp.x) * 0.5f, (hp.y - bp.y) * 0.5f);

            const double nb_dir[2] = { np.x - bp.x, np.y - bp.y };
            const double pb_dir[2] = { pp.x - bp.x, pp.y - bp.y };

            double max_radius_len = np_len > pp_len ? pp_len : np_len;
            max_radius_len = max_radius_len > hp_len ? hp_len : max_radius_len; 

            double actual_radius = (max_radius_len / 100.0) * radius;

            double ax0, ay0, ax1, ay1;
            VectorExtend(
                bp.x, bp.y,
                np.x - bp.x,
                np.y - bp.y,
                actual_radius,
                actual_radius,
                &ax0, &ay0
            );

            VectorExtend(
                bp.x, bp.y,
                pp.x - bp.x,
                pp.y - bp.y,
                actual_radius,
                actual_radius,
                &ax1, &ay1
            );

            double rad0 = atan2(ay0 - bp.y, ax0 - bp.x);
            double rad1 = atan2(ay1 - bp.y, ax1 - bp.x);

            double radmin = rad0 < rad1 ? rad0 : rad1;
            double radmax = rad0 > rad1 ? rad0 : rad1;

            if (radmax - radmin > M_PI)
            {
                radmin = radmin < 0 ? radmin + M_PI * 2.0 : radmin;
                radmax = radmax < 0 ? radmax + M_PI * 2.0 : radmax;
            }

            if (filled_angles)
                draw_list->PathLineTo(bp);

            draw_list->PathArcTo(bp, (float)actual_radius, (float)radmin, (float)radmax);

            if (filled_angles)
                draw_list->PathFillConvex(0xffffffff);
            else
                draw_list->PathStroke(0xffffffff);

            // points
            char buf[3];
            buf[0] = 'A' + i;
            buf[1] = '\0';

            ImVec2 text_size = ImGui::CalcTextSize(buf);

            double rpx, rpy;
            VectorExtend(
                bp.x, bp.y,
                -(nb_dir[0] + pb_dir[0]), -(nb_dir[1] + pb_dir[1]),
                text_size.x * 1.5f, text_size.y * 1.5f,
                &rpx, &rpy
            );

            draw_list->AddText({ (float)rpx - text_size.x * 0.5f, (float)rpy - text_size.y * 0.5f }, ImGui::GetColorU32(ImGuiCol_Text), buf);

            // median letters
            if (show_medians)
            {
                const ImVec2 mp = ImPlot::PlotToPixels({ data.medians_x[i], data.medians_y[i] });
                VectorExtend(
                    mp.x, mp.y,
                    mp.x - bp.x, mp.y - bp.y,
                    text_size.x * 2.f, text_size.y * 2.f,
                    &rpx, &rpy
                );

                buf[0] = 'M';
                buf[1] = '1' + i;
                buf[2] = '\0';
                draw_list->AddText({ (float)rpx - text_size.x * 0.5f, (float)rpy - text_size.y * 0.5f }, ImGui::GetColorU32(ImGuiCol_Text), buf);
            }

            // height letters
            if (show_heights)
            {
                VectorExtend(
                    hp.x, hp.y,
                    hp.x - bp.x, hp.y - bp.y,
                    text_size.x * 2.f, text_size.y * 2.f,
                    &rpx, &rpy
                );

                buf[0] = 'H';
                buf[1] = '1' + i;
                buf[2] = '\0';
                draw_list->AddText({ (float)rpx - text_size.x * 0.5f, (float)rpy - text_size.y * 0.5f }, ImGui::GetColorU32(ImGuiCol_Text), buf);
            }
        }   

        ImPlot::PopPlotClipRect();

        ImPlot::PlotLine("##TriangleSides", data.plot_x, data.plot_y, 3, ImPlotLineFlags_Loop);

        for (int i = 0; i < 3; i++)
        {
            if (ImPlot::DragPoint(i, &mypoints[i][0], &mypoints[i][1], ImGui::GetStyleColorVec4(ImGuiCol_Text)))
                wants_data_update = true;
        }

        if (show_medians)
        {
            ImPlot::PlotScatter("##MedianPoints", data.medians_x, data.medians_y, 3);

            for (int i = 0; i < 3; i++)
            {
                double tmp_x[2];
                double tmp_y[2];

                tmp_x[0] = data.plot_x[i];
                tmp_x[1] = data.medians_x[i];

                tmp_y[0] = data.plot_y[i];
                tmp_y[1] = data.medians_y[i];

                char buf[32];
                snprintf(buf, sizeof(buf), "##Median %c", 'A' + i);
                ImPlot::PlotLine(buf, tmp_x, tmp_y, 2);
            }

            ImPlot::PlotScatter("##Centroid", &data.centroid_x, &data.centroid_y, 1);
        }

        if (show_heights)
        {
            for (int i = 0; i < 3; i++)
            {
                double tmp_x[2];
                double tmp_y[2];

                tmp_x[0] = data.plot_x[i];
                tmp_x[1] = data.heights_x[i];

                tmp_y[0] = data.plot_y[i];
                tmp_y[1] = data.heights_y[i];

                char buf[32];
                snprintf(buf, sizeof(buf), "##Height %c", 'A' + i);
                ImPlot::PlotLine(buf, tmp_x, tmp_y, 2);
            }

            ImPlot::PlotScatter("##HeightBases", data.heights_x, data.heights_y, 3);
        }

        /*ImPlot::PlotText("A", data.plot_x[0], data.plot_y[0]);
        ImPlot::PlotText("B", data.plot_x[1], data.plot_y[1]);
        ImPlot::PlotText("C", data.plot_x[2], data.plot_y[2]);*/

        char buf[16];

        // lengths
        if (show_lengths)
        {
            for (int i = 0; i < 3; i++)
            {
                const int ni = i == 2 ? 0 : i + 1;

                snprintf(buf, sizeof(buf), "%.3f", data.lengths[i]);
                ImPlot::PlotText(buf,
                    (data.plot_x[i] + data.plot_x[ni]) * 0.5,
                    (data.plot_y[i] + data.plot_y[ni]) * 0.5
                );
            }
        }
        ImPlot::EndPlot();
    }

    ImGui::TextDisabled("Points:");
    for (int i = 0; i < 3; i++)
        ImGui::Text("%c = (%g, %g)", 'A' + i, data.plot_x[i], data.plot_y[i]);
    ImGui::Spacing();

    ImGui::TextDisabled("Area:");
    ImGui::Text("S = %g", data.area);

    ImGui::TextDisabled("Angles:");
    for (int i = 0; i < 3; i++)
        ImGui::Text("%c = %g rad (%g deg)", 'A' + i, data.angles[i], data.angles[i] * (180.0 / M_PI));
    ImGui::Spacing();

    ImGui::TextDisabled("Side Equations:");
    for (int i = 0; i < 3; i++)
    {
        char buf[64];
        FormatEquation(buf, sizeof(buf), data.side_equations[i]);
        ImGui::TextUnformatted(buf);
    }
    ImGui::Spacing();

    ImGui::TextDisabled("Lengths:");
    ImGui::Text("|AB| = %g", data.lengths[0]);
    ImGui::Text("|BC| = %g", data.lengths[1]);
    ImGui::Text("|AC| = %g", data.lengths[2]);
    ImGui::Spacing();

    ImGui::TextDisabled("Height Bases:");
    for (int i = 0; i < 3; i++)
        ImGui::Text("%cH%i = (%g, %g)", 'A' + i, i, data.heights_x[i], data.heights_y[i]);
    ImGui::Spacing();

    ImGui::TextDisabled("Height Lengths:");
    for (int i = 0; i < 3; i++)
        ImGui::Text("|%cH%i| = %g", 'A' + i, i, data.heights[i]);
    ImGui::Spacing();

    ImGui::TextDisabled("Height Equations:");
    for (int i = 0; i < 3; i++)
    {
        char buf[64];
        FormatEquation(buf, sizeof(buf), data.height_equations[i]);
        ImGui::TextUnformatted(buf);
    }
    ImGui::Spacing();

    ImGui::TextDisabled("Median Bases:");
    for (int i = 0; i < 3; i++)
        ImGui::Text("%cM%i = (%g, %g)", 'A' + i, i, data.medians_x[i], data.medians_y[i]);
    ImGui::Spacing();

    ImGui::TextDisabled("Median Lengths:");
    for (int i = 0; i < 3; i++)
        ImGui::Text("|%cM%i| = %g", 'A' + i, i, data.median_lengths[i]);
}

void RunMainWindow()
{
    if (ImGui::BeginTabBar("##TabBar"))
    {
        if (ImGui::BeginTabItem("Triangle Inspection"))
        {
            RenderTriangleWindow();
            ImGui::EndTabItem();
        }
        
        if (ImGui::BeginTabItem("Calculations"))
        {
            ImGui::TextUnformatted("Not implemented yet!");
            ImGui::EndTabItem();
        }

        ImGui::EndTabBar();
    }
}