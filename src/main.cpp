#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;
  auto start = chrono::system_clock::now();
  h.onMessage([&mpc,&start](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          
          v = v*0.44704; //convert speed unit to m/s

          /*
          * TODO: Calculate steering angle and throttle using MP0.44704C.
          *
          * Both are in between [-1, 1].
          *
          */
          double steer_value;
          double throttle_value;
// ============================= mycode ===========================

          // the car's state after the latency 
          double latency = 0.15;// computation time + latency
          px = px+v*cos(psi)*latency;
          py = py+v*sin(psi)*latency;
          psi = psi+v*mpc.last_delta/2.67*latency;
          v = v+mpc.last_a*latency;

          // convert the sensed points to vehicle's coordinate
          int num_pts = ptsx.size();
          for(int i=0; i<num_pts; i++){
            double car_x = (ptsx[i]-px)*cos(psi)+(ptsy[i]-py)*sin(psi);
            double car_y = (ptsy[i]-py)*cos(psi)-(ptsx[i]-px)*sin(psi);
            ptsx[i] = car_x;
            ptsy[i] = car_y;
          }

          // get the coefficent
          Eigen::VectorXd x_pts(num_pts);
          Eigen::VectorXd y_pts(num_pts);
          for(int i=0; i<num_pts;i++){
            x_pts[i] = ptsx[i];
            y_pts[i] = ptsy[i];
          }
          Eigen::VectorXd coeffs = polyfit(x_pts, y_pts, 3);

          // get the cte and epsi
          double cte = 0 - polyeval(coeffs,0);
          double epsi = 0 - atan(coeffs[1]);

          // set state 
          Eigen::VectorXd state(6);
          state << 0, 0, 0, v, cte, epsi;
          
          //get mpc values
          auto mpc_vars = mpc.Solve(state, coeffs);

          // set steer and throttle
          steer_value = -mpc_vars[0]/deg2rad(25); // convert to [-1, 1]
          throttle_value = mpc_vars[1]/0.44704; // convert to [-1, 1]
          
// ================================================================
          
          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          size_t num_mpc_pts = (mpc_vars.size())/2-1;
          for(size_t i=0; i< num_mpc_pts; i++){
            mpc_x_vals.push_back(mpc_vars[2+i]);
            mpc_y_vals.push_back(mpc_vars[2+i+num_mpc_pts]);
          }

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;
          
          //Display the waypoints/reference line
          // vector<double> next_x_vals(ptsx.begin(), ptsx.end());
          // vector<double> next_y_vals(ptsy.begin(), ptsy.end());

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = ptsx;
          msgJson["next_y"] = ptsy;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl<<std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

          // get the computation time
          auto end = chrono::system_clock::now();
          auto dur = chrono::duration_cast<chrono::milliseconds>(end - start);
          std::cout<<"computation time: "<< dur.count() <<"  ms"<<std::endl;
          start = end;
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
