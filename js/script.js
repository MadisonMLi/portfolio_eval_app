function handle_submit() {

    var request = new XMLHttpRequest();
    request.open("POST", "http://localhost:5115/" + button_european.value);
    request.setRequestHeader("Content-Type", "application/json");

    var price_output = document.getElementById("price_output");
    var se_output = document.getElementById("se_output");
    var delta_output = document.getElementById("delta_output");
    var gamma_output = document.getElementById("gamma_output");
    var vega_output = document.getElementById("vega_output");
    var theta_output = document.getElementById("theta_output");
    var rho_output = document.getElementById("rho_output");

    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            var output = JSON.parse(request.responseText);
            console.log(output);

            price_output.innerHTML = output.price;
            se_output.innerHTML = output.se;
            delta_output.innerHTML = output.delta;
            gamma_output.innerHTML = output.gamma;
            vega_output.innerHTML = output.vega;
            theta_output.innerHTML = output.theta;
            rho_output.innerHTML = output.rho;
        }
    };

    var form_data = new FormData(document.forms.form_european);

    var object = {};
    form_data.forEach(function(value, key) {
        if (value === 'True') {
            value = true;
        } 

        if (value === 'False') {
            value = false;
        }

        object[key] = value;
    });
    var json = JSON.stringify(object);

    console.log(json)
    request.send(json);

}

const button_european = document.getElementById("european_button");
button_european.addEventListener("click", handle_submit);

var radios = document.forms["option_select"].elements["option_name"];
for(var i = 0, max = radios.length; i < max; i++) {
    radios[i].onclick = function() {
        var option_display = document.getElementById("option_select_display");
        button_european.value = this.value;
        option_display.innerHTML = button_european.value + " Option";
        button_european.innerHTML = "Simulate " + this.value;
        console.log(button_european.value);

        var IsCall_button = document.getElementById("IsCall_input");
        if (button_european.value == "Range" & IsCall_button.style.display == "block") {
            IsCall_button.style.display = "none";
        }
        else {
            IsCall_button.style.display = "block";
        }

        var strike_button = document.getElementById("strike_input");
        if (button_european.value == "Range" & strike_button.style.display == "block") {
            strike_button.style.display = "none";
        }
        else {
            strike_button.style.display = "block";
        }

        var payout_button = document.getElementById("payout_input");
        if (button_european.value == "Digital" & payout_button.style.display == "none") {
            payout_button.style.display = "block";
        }
        else {
            payout_button.style.display = "none";
        }

        var barrierlevel_button = document.getElementById("barrierlevel_input");
        if (button_european.value == "Barrier" & barrierlevel_button.style.display == "none") {
            barrierlevel_button.style.display = "block";
        }
        else {
            barrierlevel_button.style.display = "none";
        }

        var knocktype_button = document.getElementById("knocktype_input");
        if (button_european.value == "Barrier" & knocktype_button.style.display == "none") {
            knocktype_button.style.display = "block";
        }
        else {
            knocktype_button.style.display = "none";
        }

    }
}


// Units class
function units_submit() {

    var request = new XMLHttpRequest();
    request.open("POST", "http://localhost:5216/" + button_european.value);
    request.setRequestHeader("Content-Type", "application/json");

    var price_output = document.getElementById("price_output");
    var se_output = document.getElementById("se_output");
    var delta_output = document.getElementById("delta_output");
    var gamma_output = document.getElementById("gamma_output");
    var vega_output = document.getElementById("vega_output");
    var theta_output = document.getElementById("theta_output");
    var rho_output = document.getElementById("rho_output");

    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            var output = JSON.parse(request.responseText);
            console.log(output);

            price_output.innerHTML = output.price;
            se_output.innerHTML = output.se;
            delta_output.innerHTML = output.delta;
            gamma_output.innerHTML = output.gamma;
            vega_output.innerHTML = output.vega;
            theta_output.innerHTML = output.theta;
            rho_output.innerHTML = output.rho;
        }
    };

    var form_data = new FormData(document.forms.form_european);

    var object = {};
    form_data.forEach(function(value, key) {
        if (value === 'True') {
            value = true;
        } 

        if (value === 'False') {
            value = false;
        }

        object[key] = value;
    });
    var json = JSON.stringify(object);

    console.log(json)
    request.send(json);

}