using Dash
using DashHtmlComponents
using DashCoreComponents

function visualize_conserved_utrs()

    app = dash()

    app.layout = html_div(id="root") do 
        html_div(id="left-column") do
            html_h1("Conserved UTRs"),
            dcc_input(id = "my-id1", value="initial value", type = "text"),
            html_div(id = "my-div1")
        end,
        html_div(id="app-container") do
            dcc_input(id = "my-id2", value="initial value", type = "text"),
            html_div(id = "my-div2")
        end
    end

    callback!(app, Output("my-div1", "children"), Input("my-id1", "value")) do input_value
        "You've entered kot 2 $(input_value)"
    end

    callback!(app, Output("my-div2", "children"), Input("my-id2", "value")) do input_value
        "You've entered kot 2 $(input_value)"
    end

    run_server(app, "0.0.0.0", 8083; debug=true)
end

visualize_conserved_utrs()
